
from ._sample import Sample

from ..utils.geom import exterior_coords
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
        tile_scale_xfactor=1.0, 
        tile_scale_yfactor=1.0,
        **kwargs
    ):

        super().__init__(*args, **kwargs)

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
        return float((self.pixels_per_bin or 0) * self.tile_scale_xfactor)

    @property
    def native_height(self):
        return float((self.pixels_per_bin or 0) * self.tile_scale_yfactor)

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

    def _get_aggregated_sample(self, factor):
        import anndata as ad
        import scipy.sparse as sp

        if factor in self._aggregated_samples:
            return self._aggregated_samples[factor]

        report(f"Aggregating bins ({self._describe_factor(factor)})...")

        cell_width = self.native_width * factor
        cell_height = self.native_height * factor

        centroids = self.gdf.geometry.centroid
        grid = pd.DataFrame(
            {
                "grid_x": np.floor(centroids.x.to_numpy() / cell_width),
                "grid_y": np.floor(centroids.y.to_numpy() / cell_height),
            },
            index=self.gdf[self.on].astype(str).to_numpy(),
        )
        grid = grid[~grid.index.duplicated()].reindex(self.adata.obs_names.astype(str))
        valid = grid.notna().all(axis=1).to_numpy()

        # Same naming as in _fetch_aggregated_data, so obs_names match the plotted cells
        cell_names = pd.Series(np.nan, index=self.adata.obs_names, dtype=object)
        cell_names[valid] = (
            "summary:"
            + grid.loc[valid, "grid_x"].astype(np.int64).astype(str)
            + ":"
            + grid.loc[valid, "grid_y"].astype(np.int64).astype(str)
            + f":{factor}"
        ).to_numpy()

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

    def _align_bbox(self, bbox, cell_width, cell_height):
        if bbox is None:
            min_x, min_y, max_x, max_y = self.gdf.total_bounds
            bbox = [min_x, max_x, min_y, max_y]

        left, right, bottom, top = bbox

        aligned_left = np.floor(left / cell_width) * cell_width
        aligned_right = np.ceil(right / cell_width) * cell_width
        aligned_bottom = np.floor(bottom / cell_height) * cell_height
        aligned_top = np.ceil(top / cell_height) * cell_height

        return [aligned_left, aligned_right, aligned_bottom, aligned_top]

    def _get_numeric_value(self, adata, gdf, colname):
        values = self.get_value(adata=adata, gdf=gdf, colname=colname)
        if not isinstance(values, pd.Series):
            values = pd.Series(values, index=gdf.index)
        return pd.to_numeric(values, errors="coerce")

    def _fetch_aggregated_data(self, bbox=None, colnames=None, factor=1):
        cell_width = self.native_width * factor
        cell_height = self.native_height * factor

        aligned_bbox = self._align_bbox(bbox, cell_width, cell_height)

        cache_key = (
            factor,
            tuple(round(value, 4) for value in aligned_bbox),
            tuple(sorted(colnames or [])),
        )
        
        if cache_key in self._multiscale_cache:
            self._multiscale_cache.move_to_end(cache_key)
            return self._multiscale_cache[cache_key].copy()

        cropped_gdf, cropped_adata = self.crop(*aligned_bbox, crop_adata=False), self.adata
        if cropped_gdf.empty:
            cols = ["x", "y", "width", "height", "angle", "xs", "ys", "obs_name"] + (list(colnames) if colnames else [])
            empty = pd.DataFrame(columns=cols)
            self._cache_result(cache_key, empty)
            return empty.copy()

        coords = self.get_coordinates(cropped_gdf)
        grid_x = np.floor(coords["x"].to_numpy() / cell_width).astype(np.int32)
        grid_y = np.floor(coords["y"].to_numpy() / cell_height).astype(np.int32)

        data = pd.DataFrame(
            {
                "grid_x": grid_x,
                "grid_y": grid_y,
                "x": (grid_x + 0.5) * cell_width,
                "y": (grid_y + 0.5) * cell_height,
            },
            index=coords.index,
        )
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

        # These cells are defined in an axis-aligned display-space grid.
        # Applying the native tissue rotation would rotate their dimensions
        # away from their center spacing, creating overlaps and gaps.
        grouped["angle"] = 0.0
        left, right = grouped.x - cell_width / 2, grouped.x + cell_width / 2
        bottom, top = grouped.y - cell_height / 2, grouped.y + cell_height / 2
        grouped["xs"] = np.column_stack((left, right, right, left, left)).tolist()
        grouped["ys"] = np.column_stack((bottom, bottom, top, top, bottom)).tolist()
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
                data["angle"] = self.native_angle
                data["xs"], data["ys"] = exterior_coords(self.gdf.geometry.loc[data.index].values)
            return data

        return self._fetch_aggregated_data(bbox=bbox, colnames=colnames, factor=factor)