
from ._sample import Sample
from .utils.viewport import BlockViewportIndex
from ..utils.progress import report

import pandas as pd
import numpy as np
from collections import OrderedDict

class Bins(Sample):

    AGGREGATED_PREFIX = 'agg:'
    AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
    MAX_CACHE_ENTRIES = 24
    CLUSTER_PREFIX = 'leiden_'

    def __init__(
        self, 
        *args, 
        pixels_per_bin=None,
        bin_size=None,
        native_angle=0.0,
        spatial_transform=None,
        tile_scale_xfactor=1.0, 
        tile_scale_yfactor=1.0,
        aggregated_prefix=None, aggregation_factors=None, max_cache_entries=None,
        cluster_prefix=None, aggregate_obs_prefix='summary:',
        **kwargs
    ):

        super().__init__(*args, **kwargs)

        for name, value in (('AGGREGATED_PREFIX', aggregated_prefix), ('AGGREGATION_FACTORS', aggregation_factors),
                            ('MAX_CACHE_ENTRIES', max_cache_entries), ('CLUSTER_PREFIX', cluster_prefix)):
            if value is not None:
                setattr(self, name, tuple(value) if name == 'AGGREGATION_FACTORS' else value)

        self.aggregate_obs_prefix = aggregate_obs_prefix
        self.spatial_transform = spatial_transform
        self._grid_cache = {}
        self._viewport_indexes = {}
        self._grid_source = self.gdf
        self._aggregation_labels = None
        self.pixels_per_bin = pixels_per_bin
        self.bin_size = bin_size or 0
        self.native_angle = float(native_angle)

        self.tile_scale_xfactor = tile_scale_xfactor
        self.tile_scale_yfactor = tile_scale_yfactor

        self._multiscale_cache = OrderedDict()

        self._current_factor = None
        self._aggregated_samples = {}
        self._aggregated_cell_names = {}
        self._clustering_factor = 1

    def crop_geom(self, left, right, top, bottom):
        half_bin_size = self.bin_size / 2
        
        if left is None: left = half_bin_size
        if right is None: right = half_bin_size
        if bottom is None: bottom = half_bin_size
        if top is None: top = half_bin_size

        if left >= right or top >= bottom:
            return self.gdf.iloc[0:0]  # empty GeoDataFrame with same schema

        return super().crop_geom(left - half_bin_size, right + half_bin_size,
                                 top - half_bin_size, bottom + half_bin_size)

    def get_unit_area_um2(self):
        return self.bin_size ** 2

    @property
    def native_width(self):
        return self.spatial_transform.bin_width if self.spatial_transform else float((self.pixels_per_bin or 0) * self.tile_scale_xfactor)

    @property
    def native_height(self):
        return self.spatial_transform.bin_height if self.spatial_transform else float((self.pixels_per_bin or 0) * self.tile_scale_yfactor)

    def get_aggregation_factor(self, resolution_level):
        if resolution_level is None:
            return 1
        if isinstance(resolution_level, str) and resolution_level.startswith(self.AGGREGATED_PREFIX):
            return max(1, int(resolution_level.split(":", 1)[1]))
        return 1

    def get_resolution_level_for_factor(self, factor):
        if factor <= 1:
            return None
        return f"{self.AGGREGATED_PREFIX}{int(factor)}"

    @staticmethod
    def _describe_factor(factor):
        return "native bins" if factor <= 1 else f"{factor}x{factor} aggregated bins"

    def _log_resolution_change(self, factor):
        if factor != self._current_factor:
            previous = "none" if self._current_factor is None else self._describe_factor(self._current_factor)
            self._current_factor = factor

    def get_forced_aggregation_factor(self, colnames):
        """Minimum display factor required by the requested cluster columns."""
        factors = []
        for colname in colnames or ():
            if not isinstance(colname, str):
                continue
            for prefix in (f"{algorithm}_" for algorithm in ("leiden", "kmeans", "dbscan", "hdbscan", "gaussian_mixture", "agglomerative")):
                if not colname.startswith(prefix):
                    continue
                suffix = colname[len(prefix):]
                if suffix.isdigit():
                    factors.append(int(suffix))
        return max(factors) if factors else None

    def aggregation_grid(self, factor):
        """One array-block definition for display, expression and selections."""
        if self._grid_source is not self.gdf:
            self._grid_cache.clear()
            self._viewport_indexes.clear()
            self._grid_source = self.gdf
        if factor not in self._grid_cache:
            if self.spatial_transform is not None:
                cols = np.floor(self.gdf['array_col'].to_numpy() / factor).astype(np.int64)
                rows = np.floor(self.gdf['array_row'].to_numpy() / factor).astype(np.int64)
                # Integer array coordinates are native centres, not boundaries.
                centres = self.spatial_transform.array_points(
                    np.column_stack(((cols + .5) * factor - .5, (rows + .5) * factor - .5)))
            else:
                # Compatibility for in-memory samples without original positions.
                coordinates = self.get_coordinates(self.gdf)
                cols = np.floor(coordinates.x.to_numpy() / (self.native_width * factor)).astype(np.int64)
                rows = np.floor(coordinates.y.to_numpy() / (self.native_height * factor)).astype(np.int64)
                centres = np.column_stack(((cols + .5) * self.native_width * factor,
                                           (rows + .5) * self.native_height * factor))
            self._grid_cache[factor] = pd.DataFrame(dict(grid_x=cols, grid_y=rows,
                x=centres[:, 0], y=centres[:, 1]), index=self.gdf.index)
        return self._grid_cache[factor]

    def aggregation_names(self, factor):
        cache = self._aggregation_labels
        if cache is None or cache[0] is not self.gdf or cache[1] is not self.adata:
            cache = self._aggregation_labels = (self.gdf, self.adata, {})
        if factor not in cache[2]:
            grid = self.aggregation_grid(factor)
            names = self.aggregate_obs_prefix + grid.grid_x.astype(str) + ':' + grid.grid_y.astype(str) + f':{factor}'
            names.index = self.gdf[self.on].astype(str)
            cache[2][factor] = names[~names.index.duplicated()].reindex(self.adata.obs_names.astype(str))
        return cache[2][factor]

    def _get_aggregated_sample(self, factor):
        import anndata as ad
        import scipy.sparse as sp

        if factor in self._aggregated_samples:
            return self._aggregated_samples[factor]

        report(f"Aggregating bins ({self._describe_factor(factor)})...")

        cell_names = self.aggregation_names(factor)
        valid = cell_names.notna().to_numpy()

        codes, uniques = pd.factorize(cell_names[valid])
        membership = sp.csr_matrix(
            (np.ones(len(codes), dtype=np.float32), (codes, np.flatnonzero(valid))),
            shape=(len(uniques), self.adata.n_obs),
        )
        aggregated_X = sp.csr_matrix(membership @ self.adata.X)

        aggregated_adata = ad.AnnData(
            X=aggregated_X,
            obs=pd.DataFrame(index=pd.Index(uniques, dtype=str)),
            var=pd.DataFrame(index=self.adata.var_names.copy()),
        )

        aggregated_sample = Sample(
            pd.DataFrame({"index": aggregated_adata.obs_names}),
            adata=aggregated_adata,
            on="index",
            **self._sample_options,
        )
        self._aggregated_samples[factor] = aggregated_sample
        self._aggregated_cell_names[factor] = cell_names
        return aggregated_sample

    def _cache_result(self, cache_key, value):
        self._multiscale_cache[cache_key] = value
        self._multiscale_cache.move_to_end(cache_key)
        while len(self._multiscale_cache) > self.MAX_CACHE_ENTRIES:
            self._multiscale_cache.popitem(last=False)

    def _get_numeric_value(self, adata, gdf, colname, raw=False):
        values = self.get_value(adata=adata, gdf=gdf, colname=colname, raw=raw)
        if not isinstance(values, pd.Series):
            values = pd.Series(values, index=gdf.index)
        return pd.to_numeric(values, errors="coerce")

    def _crop_aggregation_grid(self, factor, bbox):
        grid = self.aggregation_grid(factor)
        if bbox is None:
            return self.gdf, grid.copy()
        left, right, bottom, top = bbox
        angle = self.spatial_transform.angle if self.spatial_transform else 0.
        cell_width, cell_height = self.native_width * factor, self.native_height * factor
        half_x = (abs(np.cos(angle))*cell_width + abs(np.sin(angle))*cell_height) / 2
        half_y = (abs(np.sin(angle))*cell_width + abs(np.cos(angle))*cell_height) / 2
        index = self._viewport_indexes.get(factor)
        if index is None or index.grid is not grid:
            index = self._viewport_indexes[factor] = BlockViewportIndex(grid)
        # Expand conservatively to include floating-point boundary contacts,
        # then apply the original inequalities to the small candidate set.
        scale = max(1., abs(left), abs(right), abs(bottom), abs(top), half_x, half_y)
        margin = np.finfo(float).eps * scale * 4
        rows = index.query(left - half_x - margin, right + half_x + margin,
                           bottom - half_y - margin, top + half_y + margin)
        candidates = grid.iloc[rows]
        intersects = ((candidates.x + half_x >= left) & (candidates.x - half_x <= right)
                      & (candidates.y + half_y >= bottom) & (candidates.y - half_y <= top))
        rows = rows[intersects.to_numpy()]
        return self.gdf.iloc[rows], grid.iloc[rows].copy()

    def _fetch_aggregated_data(self, bbox=None, colnames=None, factor=1):
        cell_width = self.native_width * factor
        cell_height = self.native_height * factor

        cache_key = (factor, None if bbox is None else tuple(round(v, 4) for v in bbox),
                     tuple(sorted(colnames or [])), self.expression_scale)
        if cache_key in self._multiscale_cache:
            self._multiscale_cache.move_to_end(cache_key)
            return self._multiscale_cache[cache_key].copy()

        cropped_gdf, data = self._crop_aggregation_grid(factor, bbox)
        cropped_adata = self.adata
        if cropped_gdf.empty:
            cols = ["x", "y", "width", "height", "angle", "obs_name"] + (list(colnames) if colnames else [])
            return pd.DataFrame(columns=cols)
        grid_x, grid_y = data.grid_x.to_numpy(), data.grid_y.to_numpy()
        
        data["width"] = cell_width
        data["height"] = cell_height
        data["obs_name"] = (
            self.aggregate_obs_prefix
            + data["grid_x"].astype(str)
            + ":"
            + data["grid_y"].astype(str)
            + f":{factor}"
        )

        grouped = data.groupby(["grid_x", "grid_y"], sort=False).agg(
            x=("x", "first"),
            y=("y", "first"),
            width=("width", "first"),
            height=("height", "first"),
            obs_name=("obs_name", "first"),
        )

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
                
            for colname in colnames:
                is_categorical = colname in self.adata.obs.columns and not self.is_continous(colname)
                if colname in self.adata.var_names and self.expression_scale != "counts":
                    # Normalize each aggregate by its own total counts: summing
                    # per-bin normalized values would weight sparse bins up.
                    frame = pd.DataFrame({
                        "grid_x": grid_x, "grid_y": grid_y,
                        "counts": self._get_numeric_value(cropped_adata, cropped_gdf, colname, raw=True).values,
                        "library": self.library_values(cropped_gdf),
                    })
                    sums = frame.groupby(["grid_x", "grid_y"], sort=False)[["counts", "library"]].sum(min_count=1)
                    grouped[colname] = pd.Series(
                        self.scale_expression(sums["counts"].to_numpy(), sums["library"].to_numpy()),
                        index=sums.index)
                    continue
                if is_categorical:
                    values = self.get_value(adata=cropped_adata, gdf=cropped_gdf, colname=colname)
                else:
                    values = self._get_numeric_value(cropped_adata, cropped_gdf, colname)
                value_frame = pd.DataFrame(
                    {
                        "grid_x": grid_x,
                        "grid_y": grid_y,
                        colname: values.values,
                    },
                    index=cropped_gdf.index,
                )

                _grp = value_frame.groupby(["grid_x", "grid_y"], sort=False)[colname]
                if is_categorical:
                    grouped[colname] = _grp.agg(
                        lambda s: s.value_counts().idxmax() if len(s.value_counts()) else np.nan
                    )
                elif colname in self.adata.uns.get("mean_obs_columns", []):
                    grouped[colname] = _grp.mean()
                else:
                    grouped[colname] = _grp.sum(min_count=1)

        grouped["angle"] = self.spatial_transform.angle if self.spatial_transform else 0.0
        grouped.reset_index(drop=True, inplace=True)
        grouped.index = grouped["obs_name"]
        self._cache_result(cache_key, grouped)

        return grouped.copy()

    def fetch_data(self, bbox=None, colnames=None, return_coords=True, resolution_level=None) -> pd.DataFrame:
        factor = self.get_aggregation_factor(resolution_level)
        self._log_resolution_change(factor)
        if factor <= 1:
            data = super().fetch_data(bbox=bbox, colnames=colnames, return_coords=return_coords)
            if return_coords:
                data["width"] = self.native_width
                data["height"] = self.native_height
                data["angle"] = self.spatial_transform.angle if self.spatial_transform else 0.0
            return data

        return self._fetch_aggregated_data(bbox=bbox, colnames=colnames, factor=factor)
