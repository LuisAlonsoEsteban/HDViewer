

import numpy as np
import pandas as pd
from collections import OrderedDict
import threading

import scipy.sparse as sp
import shapely


class Sample:
    """Native observations and analysis with standalone constructor defaults.

    Applications can supply cache/preprocessing/naming options explicitly. These
    per-instance values are inherited by aggregate samples; method arguments
    override analysis defaults without depending on any application settings.
    """

    _max_y = None

    GENE_CACHE_SIZE = 8

    # How gene values are displayed; analyses (DE, clustering, groups) use counts.
    EXPRESSION_SCALES = {
        "counts": "Counts",
        "normalized": "Normalized (counts per 10k)",
        "log1p": "Log-normalized (log1p)",
    }

    def __init__(
            self, 
            gdf, 
            adata=None, 
            on='index',
            *, gene_cache_size=None, cluster_col='leiden', min_gene_cells=3,
            normalization_target_sum=1e4, marker_method='wilcoxon'
        ):
        
        self.GENE_CACHE_SIZE = self.GENE_CACHE_SIZE if gene_cache_size is None else gene_cache_size
        self.cluster_col = cluster_col
        self.min_gene_cells = min_gene_cells
        self.normalization_target_sum = normalization_target_sum
        self.marker_method = marker_method
        self._sample_options = dict(gene_cache_size=self.GENE_CACHE_SIZE, cluster_col=cluster_col,
                                    min_gene_cells=min_gene_cells, normalization_target_sum=normalization_target_sum,
                                    marker_method=marker_method)
        self.gdf = gdf
        self.adata = adata

        self.adata_processed = None

        self._step_params = {}
        self._embedding_keys = {
            "pca": "X_pca_2d",
            "umap": "X_umap",
            "tsne": "X_tsne",
        }

        self.on = on
        if self.adata is not None:
            assert self.on in self.gdf.columns, "On parameter is incorrect. It must be a column of 'gdf'"

        self.fetch_lock = threading.RLock()

        self._gene_cache = OrderedDict()
        self._gene_cache_adata = None
        self._library_sizes = (None, None)
        self.expression_scale = "counts"

        self._obs_index = None
        self._centroids = None

    @property
    def n_geometries(self):
        return self.gdf.shape[0]
    
                
    def has_adata(self):
        return self.adata is not None
    
    @property
    def expression_scale(self):
        return self._expression_scale

    @expression_scale.setter
    def expression_scale(self, scale):
        if scale not in self.EXPRESSION_SCALES:
            raise ValueError(f"Unknown expression scale {scale!r}; use one of {list(self.EXPRESSION_SCALES)}")
        self._expression_scale = scale

    def library_sizes(self, adata=None):
        """Total counts per observation of `adata` (default self.adata)."""
        if adata is not None and adata is not self.adata:
            return np.asarray(adata.X.sum(axis=1), dtype=float).ravel()
        with self.fetch_lock:
            cached_adata, sizes = self._library_sizes
            if cached_adata is not self.adata:
                sizes = np.asarray(self.adata.X.sum(axis=1), dtype=float).ravel()
                self._library_sizes = (self.adata, sizes)
            return sizes

    def scale_expression(self, counts, library_sizes):
        """Counts → the display scale. NaN (not sequenced) stays NaN; empty
        observations are 0 rather than 0/0."""
        if self.expression_scale == "counts":
            return counts
        counts = np.asarray(counts, dtype=float)
        library_sizes = np.asarray(library_sizes, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            values = counts * self.normalization_target_sum / library_sizes
        values = np.where(library_sizes == 0, 0.0, values)
        return np.log1p(values) if self.expression_scale == "log1p" else values

    def get_range(self, colname):
        if colname in self.adata.var_names:
            values = self.scale_expression(self._gene_vector(colname), self.library_sizes())
            return np.nanmin(values), np.nanmax(values)
        elif colname in self.adata.obs.columns:
            col = self.adata.obs[colname]
            if hasattr(col, 'cat'):
                col = col.cat.as_ordered()
            numeric = pd.to_numeric(col, errors='coerce')

        if numeric.notna().any():
            return numeric.min(), numeric.max()

        return 0.0, 1.0
    
    def crop_geom(self, left, right, top, bottom):
        if self.gdf.empty:
            return self.gdf.copy()

        if any(value is None for value in (left, right, top, bottom)):
            xmin, ymin, xmax, ymax = self.gdf.total_bounds
            left = xmin if left is None else left
            right = xmax if right is None else right
            bottom = ymin if bottom is None else bottom
            top = ymax if top is None else top

        # GeoPandas invalidates its cached index when geometry changes. Use the
        # same exact intersection predicate as .cx, but visit nearby candidates.
        region = shapely.box(left, bottom, right, top)
        rows = self.gdf.sindex.query(region, predicate='intersects')
        return self.gdf.iloc[np.sort(rows)].copy()

    def crop(self, left, right, top, bottom, crop_adata=False):
        cropped_gdf =  self.crop_geom(left, right, top, bottom)
        if self.has_adata() and crop_adata:
            cropped_mask = self.adata.obs_names.isin(cropped_gdf.index)
            cropped_adata = self.adata[cropped_mask, :]
            return cropped_gdf, cropped_adata
        else:
            return cropped_gdf
    
    def get_unit_area_um2(self):
        return self.gdf.geometry.area.mean()
    
    def get_coordinates(self, gdf) -> pd.DataFrame:
        if self._centroids is None or self._centroids[0] is not self.gdf:
            x, y = self._compute_centroids(self.gdf)
            self._centroids = (self.gdf, x, y)

        rows = self.gdf.index.get_indexer(gdf.index) if self.gdf.index.is_unique else None
        if rows is None or (rows < 0).any():
            x, y = self._compute_centroids(gdf)
        else:
            x, y = self._centroids[1][rows], self._centroids[2][rows]

        return pd.DataFrame({"x": x, "y": y}, index=gdf.index)

    @staticmethod
    def _compute_centroids(gdf):
        centroids = shapely.centroid(np.asarray(gdf.geometry.values, dtype=object))
        return shapely.get_x(centroids), shapely.get_y(centroids)

    def _obs_positions(self, gdf):
        cached = self._obs_index
        if cached is None or cached[0] is not self.gdf or cached[1] is not self.adata:
            obs = pd.Index(self.adata.obs_names.astype(self.gdf[self.on].dtype))
            positions = None
            if obs.is_unique and self.gdf.index.is_unique:
                positions = obs.get_indexer(self.gdf[self.on])
            names = np.asarray(self.adata.obs_names.astype(str), dtype=object)
            cached = self._obs_index = (self.gdf, self.adata, positions, names)
        positions = cached[2]
        if positions is None or gdf is self.gdf:
            return positions
        rows = self.gdf.index.get_indexer(gdf.index)
        return None if (rows < 0).any() else positions[rows]

    def _gene_vector(self, gene):
        """Values of `gene` for every obs of self.adata, as a dense vector. A column
        of the sparse matrix is slow to slice, so recent genes are cached."""
        with self.fetch_lock:
            if self._gene_cache_adata is not self.adata:
                self._gene_cache.clear()
                self._gene_cache_adata = self.adata

            if gene in self._gene_cache:
                self._gene_cache.move_to_end(gene)
                return self._gene_cache[gene]

            column = self.adata.X[:, [self.adata.var_names.get_loc(gene)]]
            values = column.toarray().ravel() if sp.issparse(column) else np.asarray(column).ravel()
            self._gene_cache[gene] = values

            while len(self._gene_cache) > self.GENE_CACHE_SIZE:
                self._gene_cache.popitem(last=False)

            return values

    def library_values(self, gdf):
        """Library size of the observation behind each row of `gdf` (NaN if none)."""
        sizes = self.library_sizes()
        positions = self._obs_positions(gdf)
        if positions is not None:
            return np.where(positions >= 0, sizes[np.maximum(positions, 0)], np.nan)
        lookup = pd.Series(sizes, index=self.adata.obs_names.astype(str))
        return self.get_obs_names(self.adata, gdf).map(lookup).to_numpy(dtype=float)

    def get_adata_vector(self, adata, colname):
        if colname in adata.var_names:
            return adata[:, colname].to_df()
        elif colname in adata.obs.columns:
            return adata.obs[colname].to_frame()
        else:
            raise ValueError(f"'{colname}' not found in adata")
         

    def get_value(self, adata=None, gdf=None, colname=None, raw=False):
        """Values of `colname` for `gdf`; genes on the display scale unless `raw`."""
        if gdf is None:
            gdf = self.gdf
        if adata is None:
            adata = self.adata

        if colname in gdf.columns:
            return gdf[colname].values.tolist()

        elif self.has_adata():
            if colname in adata.var_names and adata is self.adata:
                positions = self._obs_positions(gdf)
                if positions is not None:
                    rows = np.maximum(positions, 0)
                    values = self._gene_vector(colname)[rows]
                    if not raw:
                        values = self.scale_expression(values, self.library_sizes()[rows])
                    values = np.where(positions >= 0, values, np.nan)
                    return pd.Series(values, index=gdf[self.on].to_numpy(), name=colname)

            value_df = self.get_adata_vector(adata, colname)
            if colname in adata.var_names and not raw:
                value_df[colname] = self.scale_expression(value_df[colname].to_numpy(),
                                                          self.library_sizes(adata))

            if colname in adata.obs.columns:
                obs_series = self.get_obs_names(adata, gdf)
                result = obs_series.map(value_df[colname])
                result.index = gdf.index
                return result
            
            value_df.index.name = None  # avoid ambiguity if obs index name matches self.on
            value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
            merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')
            return merged_gdf[colname]

    def get_obs_names(self, adata=None, gdf=None):
        if gdf is None:
            gdf = self.gdf

        if adata is None:
            adata = self.adata

        if adata is None:
            return gdf[self.on].astype(str)

        if adata is self.adata:
            positions = self._obs_positions(gdf)
            if positions is not None:
                names = self._obs_index[3][np.maximum(positions, 0)]
                names = np.where(positions >= 0, names, np.nan)
                return pd.Series(names, index=gdf[self.on].to_numpy(), name="obs_name")

        obs_df = pd.DataFrame({
            "obs_name": adata.obs_names.astype(str),
            self.on: adata.obs_names.astype(gdf[self.on].dtype),
        })
        merged_gdf = gdf.reset_index().merge(obs_df, on=self.on, how='left').set_index('index')
        return merged_gdf["obs_name"]

    
    def fetch_data(self, bbox=None, colnames=None, return_coords=True, **coord_kwargs) -> pd.DataFrame:
        cropped_gdf = self.crop(*bbox, crop_adata=False) if bbox is not None else self.gdf
        cropped_adata = self.adata

        if return_coords:
            data = self.get_coordinates(cropped_gdf, **coord_kwargs)
        else:
            data = pd.DataFrame(index=cropped_gdf.index)

        data['obs_name'] = self.get_obs_names(cropped_adata, cropped_gdf).values

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                data[colname] = self.get_value(cropped_adata, cropped_gdf, colname)

        return data
    
    def get_colnames(self):
        return self.adata.var_names.tolist() + self.adata.obs.columns.tolist()
    
    def is_continous(self, colname):
        if colname in self.get_colnames():
            if colname in self.adata.var_names:
                return True
            else:
                return pd.api.types.is_numeric_dtype(self.adata.obs[colname])
        else:
            raise ValueError(f"{colname} not found in 'adata'.")

