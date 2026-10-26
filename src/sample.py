
import geopandas as gpd
import numpy as np
import pandas as pd
from collections import OrderedDict
import math
import threading

import scipy.sparse as sp
import shapely

from .utils.geom import exterior_coords
from .utils.progress import report

class Sample:

    _max_y = None
    # Genes whose values are kept as dense vectors, one float per obs each
    # (~25 MB per gene with 2 µm bins)
    GENE_CACHE_SIZE = 8

    def __init__(
            self, 
            gdf, 
            adata=None, 
            pixels_per_bin=None, 
            on='index'):
        
        self.gdf = gdf
        self.adata = adata

        self.adata_processed = None
        self._step_params = {}
        self._embedding_keys = {
            "pca": "X_pca_2d",
            "umap": "X_umap",
            "tsne": "X_tsne",
        }

        self.pixels_per_bin = pixels_per_bin

        self.on = on
        if self.adata is not None:
            assert self.on in self.gdf.columns, "On parameter is incorrect. It must be a column of 'gdf'"

        # Fetches may run in a worker thread (DynamicSourceManager.fetch_data_async):
        # this lock serialises them with each other and with cache invalidations
        self.fetch_lock = threading.RLock()
        self._gene_cache = OrderedDict()        # gene -> dense values over self.adata obs
        self._gene_cache_adata = None
        self._obs_index = None                  # (gdf, adata, obs row per gdf row, obs names)
        self._centroids = None                  # (gdf, centroid x, centroid y)

    @property
    def n_geometries(self):
        return self.gdf.shape[0]
    
    @property
    def n_not_sequenced(self):
        colname = self.adata.var_names.tolist()[0]
        return int(self.get_value(self.adata, self.gdf, colname=colname).isna().sum())
                
    def has_adata(self):
        return self.adata is not None
    
    def get_range(self, colname):
        if colname in self.adata.var_names:
            values = self._gene_vector(colname)
            return values.min(), values.max()
        elif colname in self.adata.obs.columns:
            col = self.adata.obs[colname]
            if hasattr(col, 'cat'):
                col = col.cat.as_ordered()
            numeric = pd.to_numeric(col, errors='coerce')
        if numeric.notna().any():
            return numeric.min(), numeric.max()
        return 0.0, 1.0
    
    def crop_geom(self, left, right, top, bottom):
        return self.gdf.cx[left:right, bottom:top].copy()

    def crop(self, left, right, top, bottom, crop_adata=False):
        cropped_gdf =  self.crop_geom(left, right, top, bottom)
        if crop_adata:
            cropped_mask = self.adata.obs_names.isin(cropped_gdf.index)
            cropped_adata = self.adata[cropped_mask, :]
            return cropped_gdf, cropped_adata
        else:
            return cropped_gdf
    
    def get_unit_area_um2(self):
        return self.gdf.geometry.area.mean()
    
    def get_coordinates(self, gdf) -> pd.DataFrame:
        """Centroid of each geometry (NaN for empty ones)."""
        if self._centroids is None or self._centroids[0] is not self.gdf:
            x, y = self._compute_centroids(self.gdf)
            self._centroids = (self.gdf, x, y)
        rows = self.gdf.index.get_indexer(gdf.index) if self.gdf.index.is_unique else None
        if rows is None or (rows < 0).any():           # not a subset of self.gdf
            x, y = self._compute_centroids(gdf)
        else:
            x, y = self._centroids[1][rows], self._centroids[2][rows]
        return pd.DataFrame({"x": x, "y": y}, index=gdf.index)

    @staticmethod
    def _compute_centroids(gdf):
        centroids = shapely.centroid(np.asarray(gdf.geometry.values, dtype=object))
        return shapely.get_x(centroids), shapely.get_y(centroids)

    def _obs_positions(self, gdf):
        """Row of self.adata for each row of `gdf` (-1: not sequenced), or None when
        `gdf` is not a subset of self.gdf or names are ambiguous (slow path then)."""
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

    def get_adata_vector(self, adata, colname):
        if colname in adata.var_names:
            return adata[:, colname].to_df()
        elif colname in adata.obs.columns:
            return adata.obs[colname].to_frame()
        else:
            raise ValueError(f"'{colname}' not found in adata")
         

    def get_value(self, adata=None, gdf=None, colname=None):
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
                    values = self._gene_vector(colname)[np.maximum(positions, 0)]
                    values = np.where(positions >= 0, values, np.nan)
                    return pd.Series(values, index=gdf[self.on].to_numpy(), name=colname)

            value_df = self.get_adata_vector(adata, colname)
            
            # Columnas de obs (leiden, Group 1, etc.) → map directo por obs_names
            if colname in adata.obs.columns:
                obs_series = self.get_obs_names(adata, gdf)
                result = obs_series.map(value_df[colname])
                result.index = gdf.index
                return result
            
            # Genes (adata.var_names) → merge original por self.on
            value_df.index.name = None  # avoid ambiguity if obs index name matches self.on
            value_df[self.on] = value_df.index.astype(gdf[self.on].dtype)
            merged_gdf = gdf.reset_index().merge(value_df, on=self.on, how='left').set_index('index')
            return merged_gdf[colname]
        
    def get_dummy_colname(self, dtype="numpy"):
        colname = self.get_value(colname=self.get_colnames()[0]).isna()
        if dtype == "numpy":
            return colname.to_numpy()
        elif dtype == "list":
            return colname.values.tolist()

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
        # return merged_gdf["obs_name"].astype(str)
        return merged_gdf["obs_name"]

    
    def fetch_data(self, bbox=None, colnames=None, return_coords=True, **coord_kwargs) -> pd.DataFrame:
        # The values are looked up in the whole AnnData by position: cropping it
        # first (an AnnData view) costs more than the lookup
        cropped_gdf = self.crop(*bbox, crop_adata=False) if bbox is not None else self.gdf
        cropped_adata = self.adata

        if return_coords:
            data = self.get_coordinates(cropped_gdf, **coord_kwargs)
        else:
            data = pd.DataFrame()

        data['obs_name'] = self.get_obs_names(cropped_adata, cropped_gdf).values

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                data[colname] = self.get_value(cropped_adata, cropped_gdf, colname)

        return data
    
    def get_marker_genes(self, n_genes=5):
        """Top marker genes per cluster; {} until `compute_marker_genes` has run for
        the current clustering (results of a previous clustering are not returned)."""
        if self.adata_processed is None or 'rank_genes_groups' not in self.adata_processed.uns:
            return {}
        if not self._is_cached("markers", self._step_params.get("leiden")):
            return {}
        
        result = self.adata_processed.uns['rank_genes_groups']
        clusters = result['names'].dtype.names
        markers = {}
        
        cluster_labels = self.adata_processed.obs['leiden']

        for cluster in clusters:
            genes = [result['names'][cluster][i] for i in range(n_genes)]
            scores = [result['scores'][cluster][i] for i in range(n_genes)]
            pvals = [result['pvals_adj'][cluster][i] for i in range(n_genes)]
            
            # Máscara booleana para separar las células que están dentro y fuera del cluster actual
            cluster_mask = cluster_labels == cluster
            logfcs = []
            
            for gene in genes:
                if gene in self.adata_processed.var_names:
                    # Extraemos la expresión del gen de forma segura (sea matriz densa o dispersa)
                    expr = self.adata_processed[:, gene].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    
                    # Calculamos las medias logarítmicas
                    mean_in = expr[cluster_mask].mean()
                    mean_out = expr[~cluster_mask].mean()
                    
                    # Como los datos de Scanpy ya están en log1p, aplicamos la resta directa
                    if mean_in == 0 and mean_out == 0:
                        l2fc = 0.0
                    elif mean_out == 0:
                        # Si fuera del cluster es 0, usamos pseudocuenta de seguridad para evitar división por cero
                        l2fc = np.log2((mean_in + 1e-6) / 1e-6)
                    else:
                        l2fc = mean_in - mean_out
                        
                    logfcs.append(round(float(l2fc), 2))
                else:
                    logfcs.append(0.0)
            

            # Empaquetamos los 4 elementos asegurando que todos mantengan tipos numéricos puros (float)
            markers[cluster] = list(zip(
                genes, 
                [round(float(s), 2) for s in scores], 
                [float(p) for p in pvals], 
                [float(l) for l in logfcs]
            ))
            
        return markers
    
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


    # Each processing step lists the steps whose results depend on it, so that
    # recomputing a step (e.g. with different parameters) invalidates them.
    _STEP_DEPENDENTS = {
        "pca": ("neighbors", "tsne", "pca_embedding"),
        "neighbors": ("umap", "leiden"),
        "leiden": ("markers",),
    }
    _STEP_OBSM_KEYS = {
        "pca_embedding": "X_pca_2d",
        "umap": "X_umap",
        "tsne": "X_tsne",
    }

    def _is_cached(self, step, params):
        return step in self._step_params and self._step_params[step] == params

    def _mark_computed(self, step, params):
        self._invalidate_dependents(step)
        self._step_params[step] = params

    def _invalidate_dependents(self, step):
        for dependent in self._STEP_DEPENDENTS.get(step, ()):
            self._step_params.pop(dependent, None)
            obsm_key = self._STEP_OBSM_KEYS.get(dependent)
            if obsm_key is not None and self.adata is not None:
                self.adata.obsm.pop(obsm_key, None)
            self._invalidate_dependents(dependent)

    def prepare_adata_processed(self):
        """Normalize, log-transform, select HVGs and scale. Computed only once."""
        import scanpy as sc

        if self.adata_processed is not None:
            return self.adata_processed

        report(f"Preprocessing {self.adata.n_obs} observations (normalize, HVG, scale)...")
        adata = self.adata.copy()
        adata.var_names_make_unique()

        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)

        self.adata_processed = adata
        self._step_params = {}
        return self.adata_processed

    def _ensure_pca(self, n_pcs=30):
        import scanpy as sc

        adata = self.prepare_adata_processed()
        # A PCA with at least n_pcs components can be reused by any consumer
        # asking for n_pcs (neighbors and t-SNE only use the first n_pcs).
        cached = self._step_params.get("pca")
        if cached is None or cached < n_pcs:
            report(f"Computing PCA ({n_pcs} components)...")
            sc.pp.pca(adata, n_comps=n_pcs)
            self._mark_computed("pca", n_pcs)
        return adata

    def _ensure_neighbors(self, n_pcs=30, n_neighbors=15):
        import scanpy as sc

        adata = self._ensure_pca(n_pcs)
        params = (n_pcs, n_neighbors)
        if not self._is_cached("neighbors", params):
            report(f"Computing neighbor graph (k={n_neighbors})...")
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            self._mark_computed("neighbors", params)
        return adata

    def _store_embedding(self, method, values):
        key = self._embedding_keys[method]
        full = np.full((self.adata.n_obs, 2), np.nan, dtype=float)
        obs_index = pd.Index(self.adata.obs_names.astype(str))
        processed_index = pd.Index(self.adata_processed.obs_names.astype(str))
        common_idx = obs_index.get_indexer(processed_index)
        valid = common_idx >= 0
        full[common_idx[valid], :] = values[valid, :]
        self.adata.obsm[key] = full

    def compute_embedding(self, method: str, n_pcs=30, n_neighbors=15):
        import scanpy as sc

        method = method.lower()
        if method not in self._embedding_keys:
            raise ValueError(f"Unsupported embedding method: {method}")

        key = self._embedding_keys[method]
        step = "pca_embedding" if method == "pca" else method
        params = (n_pcs, n_neighbors) if method == "umap" else n_pcs

        if self.adata_processed is not None and self._is_cached(step, params) and key in self.adata.obsm:
            report(f"Using cached {method.upper()} embedding")
            return self.adata.obsm[key]

        if method == "pca":
            adata = self._ensure_pca(n_pcs)
            values = adata.obsm["X_pca"][:, :2]
        elif method == "umap":
            adata = self._ensure_neighbors(n_pcs, n_neighbors)
            report("Running UMAP...")
            sc.tl.umap(adata)
            values = adata.obsm["X_umap"][:, :2]
        else:
            adata = self._ensure_pca(n_pcs)
            report("Running t-SNE...")
            sc.tl.tsne(adata, n_pcs=n_pcs)
            values = adata.obsm["X_tsne"][:, :2]

        self._store_embedding(method, np.asarray(values, dtype=float))
        self._mark_computed(step, params)
        return self.adata.obsm[key]

    def get_embedding_dataframe(self, method: str, resolution_level=None) -> pd.DataFrame:
        """`resolution_level` is only meaningful for Bins; it is ignored here."""
        method = method.lower()
        values = np.asarray(self.compute_embedding(method), dtype=float)
        df = pd.DataFrame(values, columns=["x", "y"], index=self.adata.obs_names.astype(str))
        df.index.name = "obs_name"
        df.reset_index(inplace=True)
        df = df.dropna(subset=["x", "y"]).copy()
        return df

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15, key_added="leiden",
                       resolution_level=None, compute_markers=True):
        """Run Leiden clustering and store the labels in `adata.obs[key_added]`.

        `resolution_level` is only meaningful for Bins; it is ignored here. With
        compute_markers=False the (slower) marker genes are left for a later
        `compute_marker_genes()`, so the clusters can be shown first.
        Returns the obs column holding the labels.
        """
        import scanpy as sc

        adata = self._ensure_neighbors(n_pcs, n_neighbors)

        if not self._is_cached("leiden", resolution):
            report(f"Running Leiden clustering (resolution={resolution})...")
            sc.tl.leiden(adata, resolution=resolution)
            self._mark_computed("leiden", resolution)
        if compute_markers:
            self.compute_marker_genes()

        common_idx = self.adata.obs_names.intersection(adata.obs_names)

        self.adata.obs[key_added] = None
        self.adata.obs.loc[common_idx, key_added] = (
            adata.obs.loc[common_idx, "leiden"]
            .astype(str)
            .values
        )
        return key_added


    def compute_marker_genes(self):
        """Rank marker genes of the current Leiden clusters (Wilcoxon); cached per clustering."""
        import scanpy as sc

        adata = self.adata_processed
        if adata is None or 'leiden' not in adata.obs:
            return
        leiden_params = self._step_params.get("leiden")
        if self._is_cached("markers", leiden_params):
            return
        report("Ranking marker genes...")
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
        self._mark_computed("markers", leiden_params)


class Bins(Sample):
    AGGREGATED_PREFIX = "agg:"
    AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
    MAX_CACHE_ENTRIES = 24
    CLUSTER_PREFIX = "leiden_"

    def __init__(self, *args, half_bin_size=None,
                 tile_scale_xfactor=1.0, tile_scale_yfactor=1.0,
                 native_angle=0.0, bin_size_um=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.half_bin_size = half_bin_size or 0
        self.native_angle = float(native_angle)
        self.bin_size_um = bin_size_um
        self.tile_scale_xfactor = tile_scale_xfactor
        self.tile_scale_yfactor = tile_scale_yfactor
        self._multiscale_cache = OrderedDict()

        self._current_factor = None
        # Per aggregation factor: a Sample holding the aggregated AnnData (so its
        # PCA/neighbors are cached across clustering runs) and, for every native
        # obs, the name of the aggregated cell it belongs to.
        self._aggregated_samples = {}
        self._aggregated_cell_names = {}
        self._clustering_factor = 1

    def crop_geom(self, left, right, top, bottom):
        if left is None: left = self.half_bin_size
        if right is None: right = self.half_bin_size
        if bottom is None: bottom = self.half_bin_size
        if top is None: top = self.half_bin_size
        if left >= right or top >= bottom:
            return self.gdf.iloc[0:0]  # empty GeoDataFrame with same schema
        result = self.gdf.cx[left - self.half_bin_size : right + self.half_bin_size,
                              top - self.half_bin_size : bottom + self.half_bin_size]
        return result.copy()

    def get_unit_area_um2(self):
        return (self.half_bin_size * 2) ** 2

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
        """Aggregation factor a 'leiden_{factor}' column was computed at, or None."""
        if not isinstance(colname, str) or not colname.startswith(self.CLUSTER_PREFIX):
            return None
        suffix = colname[len(self.CLUSTER_PREFIX):]
        if not suffix.isdigit() or self.adata is None or colname not in self.adata.obs.columns:
            return None
        return int(suffix)

    def get_forced_aggregation_factor(self, colnames):
        """Factor the view must be pinned to while showing a clustering column, or None."""
        for colname in colnames or ():
            factor = self.get_cluster_factor(colname)
            if factor is not None:
                return factor
        return None

    def _get_aggregated_sample(self, factor):
        """Sum raw counts of the native bins falling in each grid cell of size `factor`."""
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

    def get_embedding_dataframe(self, method: str, resolution_level=None) -> pd.DataFrame:
        """Embed the obs at `resolution_level`.

        For aggregated levels the embedding is computed on the aggregated AnnData
        (sharing its PCA/neighbors with clustering), and obs_names are the
        'summary:...' names of the plotted aggregated cells.
        """
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
                # Render the same footprints used for spatial operations. A
                # single angle/width/height cannot express anisotropic scaling
                # of a rotated grid (the transformed bins are parallelograms).
                data["xs"], data["ys"] = exterior_coords(self.gdf.geometry.loc[data.index].values)
            return data

        return self._fetch_aggregated_data(bbox=bbox, colnames=colnames, factor=factor)
    
class Spots(Sample):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_unit_area_um2(self):
        return 3.14 * (55.0 / 2)**2



class Cells(Sample):
    # Levels of detail: level k draws the polygons simplified with tolerance
    # LOD_BASE_TOLERANCE * 2**k (figure units); None draws the original polygons
    LOD_PREFIX = "lod:"
    LOD_BASE_TOLERANCE = 1.0
    LOD_MAX_LEVEL = 8
    # Largest outline deviation allowed, in screen pixels
    LOD_TOLERANCE_PX = 1.0
    LOD_HYSTERESIS = 0.15

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Tolerance -> the polygons of self.gdf simplified with it, in the same order.
        # Simplifying the whole slide takes ~1 s, so each level is computed once.
        self._simplified = {}

    def get_lod_level(self, pixel_size, current=None):
        """Resolution level for a view where one screen pixel spans `pixel_size`
        figure units. Near a boundary `current` is kept, so that small zooms do not
        reload the layer back and forth."""
        tolerance = pixel_size * self.LOD_TOLERANCE_PX
        if not math.isfinite(tolerance) or tolerance <= 0:
            return current
        # Continuous level: k when tolerance == LOD_BASE_TOLERANCE * 2**k (native: -1)
        value = math.log2(tolerance / self.LOD_BASE_TOLERANCE)
        k = self.get_lod_index(current)
        if not (k - self.LOD_HYSTERESIS <= value < k + 1 + self.LOD_HYSTERESIS):
            k = max(-1, min(math.floor(value), self.LOD_MAX_LEVEL))
        return None if k < 0 else f"{self.LOD_PREFIX}{k}"

    def get_lod_index(self, resolution_level):
        if isinstance(resolution_level, str) and resolution_level.startswith(self.LOD_PREFIX):
            return int(resolution_level[len(self.LOD_PREFIX):])
        return -1

    def get_tolerance(self, resolution_level):
        k = self.get_lod_index(resolution_level)
        return 0.0 if k < 0 else self.LOD_BASE_TOLERANCE * 2 ** k

    def get_simplified(self, tolerance):
        """Every cell simplified with `tolerance`, aligned with self.gdf."""
        if tolerance not in self._simplified:
            geoms = np.asarray(self.gdf.geometry.values, dtype=object)
            # preserve_topology keeps every cell a valid, non-empty polygon
            self._simplified[tolerance] = shapely.simplify(geoms, tolerance, preserve_topology=True)
        return self._simplified[tolerance]

    def get_coordinates(self, gdf, tolerance=0.0) -> pd.DataFrame:
        """Exterior ring vertices of each cell (patches), simplified with `tolerance`."""
        if tolerance > 0:
            geoms = self.get_simplified(tolerance)[self.gdf.index.get_indexer(gdf.index)]
        else:
            geoms = np.asarray(gdf.geometry.values, dtype=object)
        xs, ys = exterior_coords(geoms)
        return pd.DataFrame({"x": xs, "y": ys}, index=gdf.index)

    def fetch_data(self, bbox=None, colnames=None, return_coords=True, resolution_level=None) -> pd.DataFrame:
        return super().fetch_data(bbox=bbox, colnames=colnames, return_coords=return_coords,
                                  tolerance=self.get_tolerance(resolution_level))
