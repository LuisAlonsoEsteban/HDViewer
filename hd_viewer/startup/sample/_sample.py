

import numpy as np
import pandas as pd
from collections import OrderedDict
import threading

import scipy.sparse as sp
import shapely

from ..utils.progress import report


class Sample:
    """Native observations and analysis with standalone constructor defaults.

    Applications can supply cache/preprocessing/naming options explicitly. These
    per-instance values are inherited by aggregate samples; method arguments
    override analysis defaults without depending on any application settings.
    """

    _max_y = None

    GENE_CACHE_SIZE = 8

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

        self._obs_index = None
        self._centroids = None

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
            data = pd.DataFrame()

        data['obs_name'] = self.get_obs_names(cropped_adata, cropped_gdf).values

        if colnames is not None:
            if isinstance(colnames, str):
                colnames = [colnames]
            for colname in colnames:
                data[colname] = self.get_value(cropped_adata, cropped_gdf, colname)

        return data
    
    def get_marker_genes(self, n_genes=5):
        if self.adata_processed is None or 'rank_genes_groups' not in self.adata_processed.uns:
            return {}
        if not self._is_cached("markers", self._step_params.get("leiden")):
            return {}
        
        result = self.adata_processed.uns['rank_genes_groups']
        clusters = result['names'].dtype.names
        markers = {}
        
        cluster_labels = self.adata_processed.obs[self.cluster_col]

        for cluster in clusters:
            count = min(n_genes, len(result["names"]))
            genes = [result['names'][cluster][i] for i in range(count)]
            scores = [result['scores'][cluster][i] for i in range(count)]
            pvals = [result['pvals_adj'][cluster][i] for i in range(count)]
            
            cluster_mask = cluster_labels == cluster
            logfcs = []
            
            for gene in genes:
                if gene in self.adata_processed.var_names:
                    expr = self.adata_processed[:, gene].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    
                    mean_in = expr[cluster_mask].mean()
                    mean_out = expr[~cluster_mask].mean()
                    
                    if mean_in == 0 and mean_out == 0:
                        l2fc = 0.0
                    elif mean_out == 0:
                        l2fc = np.log2((mean_in + 1e-6) / 1e-6)
                    else:
                        l2fc = mean_in - mean_out
                        
                    logfcs.append(round(float(l2fc), 2))
                else:
                    logfcs.append(0.0)
            

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
            if obsm_key is not None and self.has_adata():
                self.adata.obsm.pop(obsm_key, None)
            self._invalidate_dependents(dependent)

    def prepare_adata_processed(self):
        import scanpy as sc

        if self.adata_processed is not None:
            return self.adata_processed

        report(f"Preprocessing {self.adata.n_obs} observations (normalize, HVG, scale)...")
        adata = self.adata.copy()
        adata.var_names_make_unique()

        sc.pp.filter_genes(adata, min_cells=self.min_gene_cells)
        sc.pp.normalize_total(adata, target_sum=self.normalization_target_sum)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable].copy()

        self.adata_processed = adata
        self._step_params = {}
        return self.adata_processed

    def _ensure_pca(self, n_pcs=30):
        import scanpy as sc

        adata = self.prepare_adata_processed()
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

    def get_embedding_dataframe(self, method: str, resolution_level=None, *, n_pcs=30, n_neighbors=15) -> pd.DataFrame:
        method = method.lower()
        values = np.asarray(self.compute_embedding(method, n_pcs=n_pcs, n_neighbors=n_neighbors), dtype=float)
        df = pd.DataFrame(values, columns=["x", "y"], index=self.adata.obs_names.astype(str))
        df.index.name = "obs_name"
        df.reset_index(inplace=True)
        df = df.dropna(subset=["x", "y"]).copy()
        return df

    def run_clustering(
        self, 
        resolution=0.9,
        n_pcs=30,
        n_neighbors=15,
        key_added=None,
        compute_markers=False,
        *, algorithm="leiden", n_clusters=8, eps=0.5, min_samples=5,
        min_cluster_size=5, covariance_type="full", linkage="ward", random_state=0,
        **kwargs
    ):
        """Run expression clustering and store the labels in `adata.obs[key_added]`.

        `resolution_level` is only meaningful for Bins; it is ignored here. With
        compute_markers=False the (slower) marker genes are left for a later
        `compute_marker_genes()`, so the clusters can be shown first.
        Returns the obs column holding the labels.
        """
        import scanpy as sc

        key_added = (self.cluster_col if algorithm == "leiden" else algorithm) if key_added is None else key_added
        
        from ...board.operations.layer.clustering import ALGORITHMS, fit_clusters
        
        if algorithm not in ALGORITHMS:
            raise ValueError(f"Unsupported clustering algorithm: {algorithm}")

        if algorithm == "leiden":
            adata = self._ensure_neighbors(n_pcs, n_neighbors)
            params = (algorithm, resolution, n_pcs, n_neighbors, random_state)
        else:
            adata = self._ensure_pca(n_pcs)
            
            # Only parameters relevant to the selected algorithm affect its cache.
            specific = {
                "kmeans": (n_clusters, random_state),
                "dbscan": (eps, min_samples),
                "hdbscan": (min_cluster_size, min_samples),
                "gaussian_mixture": (n_clusters, covariance_type, random_state),
                "agglomerative": (n_clusters, linkage),
            }
            params = (algorithm, n_pcs, specific[algorithm])

        if not self._is_cached("leiden", params):
            report(f"Running {algorithm} clustering...")
            if algorithm == "leiden":
                sc.tl.leiden(adata, resolution=resolution, key_added=self.cluster_col,
                             random_state=random_state)
            else:
                labels = fit_clusters(np.asarray(adata.obsm["X_pca"])[:, :n_pcs], algorithm,
                    n_clusters=n_clusters, eps=eps, min_samples=min_samples,
                    min_cluster_size=min_cluster_size, covariance_type=covariance_type,
                    linkage=linkage, random_state=random_state)

                adata.obs[self.cluster_col] = pd.Categorical(
                    ["Noise" if label < 0 else str(label) for label in labels])

            self._mark_computed("leiden", params)
            adata.uns.pop("rank_genes_groups", None)
        
        self.adata.uns["clustering"] = {"algorithm": algorithm, "key": key_added}
        if compute_markers:
            self.compute_marker_genes()

        common_idx = self.adata.obs_names.intersection(adata.obs_names)

        self.adata.obs[key_added] = None
        self.adata.obs.loc[common_idx, key_added] = (
            adata.obs.loc[common_idx, self.cluster_col]
            .astype(str)
            .values
        )
        return key_added


    def compute_marker_genes(self):
        import scanpy as sc

        adata = self.adata_processed
        if adata is None or self.cluster_col not in adata.obs:
            return

        leiden_params = self._step_params.get("leiden")
        if self._is_cached("markers", leiden_params):
            return

        report("Ranking marker genes...")

        counts = adata.obs[self.cluster_col].value_counts()
        groups = [str(label) for label, count in counts.items()
                  if str(label) != "Noise" and count >= 2]
        adata.uns.pop("rank_genes_groups", None)

        if groups and len(counts[counts > 0]) > 1:
            sc.tl.rank_genes_groups(adata, groupby=self.cluster_col, method=self.marker_method,
                                   groups=groups)
        
        self._mark_computed("markers", leiden_params)

