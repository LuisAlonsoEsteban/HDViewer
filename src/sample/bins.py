
from ._sample import Sample

from ..utils.progress import report

import pandas as pd
import numpy as np
from collections import OrderedDict

class Bins(Sample):

    AGGREGATED_PREFIX = "agg:"
    AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
    MAX_CACHE_ENTRIES = 24
    CLUSTER_PREFIX = "leiden_"

    def __init__(
        self, 
        *args, 
        pixels_per_bin=None,
        bin_size=None,
        native_angle=0.0,
        spatial_transform=None,
        tile_scale_xfactor=1.0, 
        tile_scale_yfactor=1.0,
        **kwargs
    ):

        super().__init__(*args, **kwargs)

        self.spatial_transform = spatial_transform
        self._grid_cache = {}
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

        result = self.gdf.cx[left - half_bin_size : right + half_bin_size,
                              top - half_bin_size : bottom + half_bin_size]

        return result.copy()

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

    def get_cluster_colname(self, factor):
        return f"{self.CLUSTER_PREFIX}{int(factor)}"

    def get_cluster_factor(self, colname):
        if not isinstance(colname, str) or not colname.startswith(self.CLUSTER_PREFIX):
            return None
        suffix = colname[len(self.CLUSTER_PREFIX):]
        if not suffix.isdigit() or self.adata is None or colname not in self.adata.obs.columns:
            return None
        return int(suffix)

    def get_forced_aggregation_factor(self, colnames):
        for colname in colnames or ():
            factor = self.get_cluster_factor(colname)
            if factor is not None:
                return factor
        return None

    def aggregation_grid(self, factor):
        """One array-block definition for display, expression and selections."""
        if self._grid_source is not self.gdf:
            self._grid_cache.clear()
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
            names = 'summary:' + grid.grid_x.astype(str) + ':' + grid.grid_y.astype(str) + f':{factor}'
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
        )
        self._aggregated_samples[factor] = aggregated_sample
        self._aggregated_cell_names[factor] = cell_names
        return aggregated_sample

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15, key_added=None, resolution_level=None,
                       compute_markers=True):
        """Cluster at `resolution_level` and store labels in `adata.obs['leiden_{factor}']`.

        For aggregated levels every native bin gets the label of the aggregated
        cell containing it, so the column can be looked up like any other obs column.
        """
        factor = self.get_aggregation_factor(resolution_level)
        key_added = key_added or self.get_cluster_colname(factor)

        if factor <= 1:
            super().run_clustering(resolution=resolution, n_pcs=n_pcs, n_neighbors=n_neighbors, key_added=key_added,
                                   compute_markers=compute_markers)
        else:
            aggregated_sample = self._get_aggregated_sample(factor)
            aggregated_sample.run_clustering(resolution=resolution, n_pcs=n_pcs, n_neighbors=n_neighbors,
                                             compute_markers=compute_markers)
            labels = aggregated_sample.adata.obs["leiden"].astype(str)
            mapped = self._aggregated_cell_names[factor].map(labels)
            self.adata.obs[key_added] = mapped.where(mapped.notna(), None).astype(object)

        self._clustering_factor = factor
        # Drop cached aggregated tiles holding stale labels of a previous run
        with self.fetch_lock:
            for cache_key in [k for k in self._multiscale_cache if key_added in k[2]]:
                del self._multiscale_cache[cache_key]
        return key_added

    def get_embedding_dataframe(self, method: str, resolution_level=None):
        factor = self.get_aggregation_factor(resolution_level)
        if factor <= 1:
            return super().get_embedding_dataframe(method)
        return self._get_aggregated_sample(factor).get_embedding_dataframe(method)

    def compute_marker_genes(self):
        if self._clustering_factor > 1:
            return self._aggregated_samples[self._clustering_factor].compute_marker_genes()
        return super().compute_marker_genes()

    def get_marker_genes(self, n_genes=5):
        if self._clustering_factor > 1:
            return self._aggregated_samples[self._clustering_factor].get_marker_genes(n_genes=n_genes)
        return super().get_marker_genes(n_genes=n_genes)

    def _cache_result(self, cache_key, value):
        self._multiscale_cache[cache_key] = value
        self._multiscale_cache.move_to_end(cache_key)
        while len(self._multiscale_cache) > self.MAX_CACHE_ENTRIES:
            self._multiscale_cache.popitem(last=False)

    def _get_numeric_value(self, adata, gdf, colname):
        values = self.get_value(adata=adata, gdf=gdf, colname=colname)
        if not isinstance(values, pd.Series):
            values = pd.Series(values, index=gdf.index)
        return pd.to_numeric(values, errors="coerce")

    def _fetch_aggregated_data(self, bbox=None, colnames=None, factor=1):
        cell_width = self.native_width * factor
        cell_height = self.native_height * factor

        cache_key = (factor, None if bbox is None else tuple(round(v, 4) for v in bbox),
                     tuple(sorted(colnames or [])))
        if cache_key in self._multiscale_cache:
            self._multiscale_cache.move_to_end(cache_key)
            return self._multiscale_cache[cache_key].copy()

        grid = self.aggregation_grid(factor)
        if bbox is None:
            mask = np.ones(len(grid), dtype=bool)
        else:
            left, right, bottom, top = bbox
            angle = self.spatial_transform.angle if self.spatial_transform else 0.
            half_x = (abs(np.cos(angle))*cell_width + abs(np.sin(angle))*cell_height) / 2
            half_y = (abs(np.sin(angle))*cell_width + abs(np.cos(angle))*cell_height) / 2
            mask = ((grid.x + half_x >= left) & (grid.x - half_x <= right)
                    & (grid.y + half_y >= bottom) & (grid.y - half_y <= top))
        cropped_gdf, cropped_adata = self.gdf.loc[mask], self.adata
        if cropped_gdf.empty:
            cols = ["x", "y", "width", "height", "angle", "obs_name"] + (list(colnames) if colnames else [])
            return pd.DataFrame(columns=cols)
        data = grid.loc[mask].copy()
        grid_x, grid_y = data.grid_x.to_numpy(), data.grid_y.to_numpy()
        
        data["width"] = cell_width
        data["height"] = cell_height
        data["obs_name"] = (
            "summary:"
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
