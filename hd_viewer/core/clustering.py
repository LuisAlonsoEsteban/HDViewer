"""Clustering of a sample's observations, without any UI.

`ClusteringState` runs the algorithms (on native or aggregated bins), keeps the
cluster columns, their ordered labels and assigned colours, which clusters are
shown, and the marker genes. Its PCA / neighbor caches live on the sample, so
embeddings reuse them. The HTTP API holds one per layer view.
"""
import numpy as np
import pandas as pd

from .. import settings
from ..settings.config import clustering_kwargs
from ..startup.utils.progress import report

CLUSTER_PALETTE = settings.CLUSTER_PALETTE

SUPPORTED_ALGORITHMS = (
    "leiden",
    "kmeans",
    "dbscan",
    "hdbscan",
    "gaussian_mixture",
    "agglomerative",
)


class ClusteringState:
    def __init__(self, sample=None):
        self._sample = sample
        self._cluster_columns: dict[str, str] = {}
        self._active_cluster_col: str | None = None
        self._cluster_labels: list[str] = []   # ordered unique cluster ids for the active result
        self._color_map: dict[str, str] = {}  # cluster_id → hex color
        self._active_clusters: set[str] = set()
        self._bound_sample = None
        # Cluster column whose marker genes are ranked (None: none are current)
        self._markers_col: str | None = None

    @property
    def _base_sample(self):
        """The layer's sample (aggregated samples are derived from it)."""
        return self._sample

    @property
    def sample(self):
        return self._bound_sample if self._bound_sample is not None else self._base_sample

    @property
    def adata(self):
        return self.sample.adata

    @adata.setter
    def adata(self, value):
        self.sample.adata = value

    @property
    def _step_params(self) -> dict:
        """Computed steps (PCA, neighbors, embeddings, clusterings) with their
        parameters. Kept on the sample, next to the data they describe, so
        clustering and embeddings of the same sample reuse each other's PCA and
        neighbor graph, and re-clustering only recomputes the clustering."""
        return self.sample._step_params

    @_step_params.setter
    def _step_params(self, value):
        self.sample._step_params = value

    @property
    def adata_processed(self):
        return self.sample.adata_processed

    @adata_processed.setter
    def adata_processed(self, value):
        self.sample.adata_processed = value

    @property
    def min_gene_cells(self):
        return self.sample.min_gene_cells

    @property
    def normalization_target_sum(self):
        return self.sample.normalization_target_sum

    def has_adata(self):
        return self.adata is not None

    @property
    def cluster_col(self):
        return self._active_cluster_col

    @cluster_col.setter
    def cluster_col(self, value):
        self._active_cluster_col = value

    def get_cluster_col(self, algorithm: str | None = None) -> str | None:
        if algorithm is not None:
            return self._cluster_columns.get(algorithm)
        return self._active_cluster_col

    def get_cluster_colname(self, factor: int, algorithm: str | None = None) -> str:
        algorithm = algorithm or "leiden"
        if factor <= 1:
            return algorithm
        return f"{algorithm}_{int(factor)}"

    def get_cluster_factor(self, colname: str, algorithm: str | None = None) -> int | None:
        if not isinstance(colname, str):
            return None
        candidates = [algorithm] if algorithm is not None else SUPPORTED_ALGORITHMS
        for candidate in candidates:
            prefix = f"{candidate}_"
            if colname.startswith(prefix):
                suffix = colname[len(prefix):]
                if suffix.isdigit():
                    return int(suffix)
        return None

    def _get_active_analysis_sample(self):
        sample = self._base_sample
        active_col = self._active_cluster_col
        if active_col is None:
            return sample
        factor = self.get_cluster_factor(active_col)
        if factor is not None and hasattr(sample, "_aggregated_samples") and factor in sample._aggregated_samples:
            return sample._aggregated_samples[factor]
        return sample

    @property
    def is_clustered(self) -> bool:
        return (
            any(column in self._base_sample.adata.obs.columns for column in self._cluster_columns.values())
            and len(self._cluster_labels) > 0
        )

    def _cluster_on_sample(self, sample, *, algorithm, resolution, n_pcs, n_neighbors,
                          key_added=None, compute_markers=False, **algorithm_options):
        """Run a single clustering algorithm on either a base sample or an aggregated one."""
        previous_sample = self._bound_sample
        self._bound_sample = sample
        try:
            key_added = algorithm if key_added is None else key_added

            if algorithm not in SUPPORTED_ALGORITHMS:
                raise ValueError(f"Unsupported clustering algorithm: {algorithm}")

            if algorithm == "leiden":
                adata = self._ensure_neighbors(n_pcs, n_neighbors)
                params = (algorithm, resolution, n_pcs, n_neighbors, algorithm_options.get("random_state", 0))
            else:
                adata = self._ensure_pca(n_pcs)
                if algorithm == "kmeans":
                    params = (algorithm, n_pcs, algorithm_options.get("n_clusters", 8), algorithm_options.get("random_state", 0))
                elif algorithm == "dbscan":
                    params = (algorithm, n_pcs, algorithm_options.get("eps", 0.5), algorithm_options.get("min_samples", 5))
                elif algorithm == "hdbscan":
                    params = (algorithm, n_pcs, algorithm_options.get("min_cluster_size", 5), algorithm_options.get("min_samples", 5))
                elif algorithm == "gaussian_mixture":
                    params = (algorithm, n_pcs, algorithm_options.get("n_clusters", 8), algorithm_options.get("covariance_type", "full"), algorithm_options.get("random_state", 0))
                elif algorithm == "agglomerative":
                    params = (algorithm, n_pcs, algorithm_options.get("n_clusters", 8), algorithm_options.get("linkage", "ward"))

            if not self._is_cached(key_added, params):
                report(f"Running {algorithm} clustering...")
                if algorithm == "leiden":
                    import scanpy as sc
                    sc.tl.leiden(adata, resolution=resolution, key_added=key_added,
                                 random_state=algorithm_options.get("random_state", 0))
                else:
                    values = np.asarray(adata.obsm["X_pca"])[:, :n_pcs]
                    if algorithm == "kmeans":
                        from sklearn.cluster import KMeans
                        estimator = KMeans(
                            n_clusters=algorithm_options.get("n_clusters", 8),
                            random_state=algorithm_options.get("random_state", 0),
                            n_init=10,
                        )
                    elif algorithm == "dbscan":
                        from sklearn.cluster import DBSCAN
                        estimator = DBSCAN(
                            eps=algorithm_options.get("eps", 0.5),
                            min_samples=algorithm_options.get("min_samples", 5),
                        )
                    elif algorithm == "hdbscan":
                        from sklearn.cluster import HDBSCAN
                        estimator = HDBSCAN(
                            min_cluster_size=algorithm_options.get("min_cluster_size", 5),
                            min_samples=algorithm_options.get("min_samples", 5),
                        )
                    elif algorithm == "gaussian_mixture":
                        from sklearn.mixture import GaussianMixture
                        estimator = GaussianMixture(
                            n_components=algorithm_options.get("n_clusters", 8),
                            covariance_type=algorithm_options.get("covariance_type", "full"),
                            random_state=algorithm_options.get("random_state", 0),
                        )
                    elif algorithm == "agglomerative":
                        from sklearn.cluster import AgglomerativeClustering
                        estimator = AgglomerativeClustering(
                            n_clusters=algorithm_options.get("n_clusters", 8),
                            linkage=algorithm_options.get("linkage", "ward"),
                        )

                    labels = estimator.fit_predict(values)
                    adata.obs[key_added] = pd.Categorical(
                        ["Noise" if label < 0 else str(label) for label in labels]
                    )

                self._mark_computed(key_added, params)
                adata.uns.pop("rank_genes_groups", None)
                # Markers ranked for the previous labels of this column are stale
                if self._markers_col == key_added:
                    self._markers_col = None

            sample.adata.uns["clustering"] = {"algorithm": algorithm, "key": key_added}
            if compute_markers:
                self.compute_markers()

            common_idx = sample.adata.obs_names.intersection(adata.obs_names)
            sample.adata.obs[key_added] = None
            sample.adata.obs.loc[common_idx, key_added] = (
                adata.obs.loc[common_idx, key_added].astype(str).values
            )
            return key_added
        finally:
            self._bound_sample = previous_sample

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15, key_added=None,
                       compute_markers=False, resolution_level=None, *, algorithm="leiden",
                       **algorithm_options):
        """Run clustering and handle bins aggregation in the manager layer."""
        sample = self._base_sample

        if hasattr(sample, "get_aggregation_factor"):
            factor = sample.get_aggregation_factor(resolution_level)
            if factor > 1:
                key_added = key_added or self.get_cluster_colname(factor, algorithm=algorithm)
                aggregated_sample = sample._get_aggregated_sample(factor)
                aggregated_key = self.get_cluster_colname(factor, algorithm=algorithm)
                self._cluster_on_sample(
                    aggregated_sample,
                    algorithm=algorithm,
                    resolution=resolution,
                    n_pcs=n_pcs,
                    n_neighbors=n_neighbors,
                    key_added=aggregated_key,
                    compute_markers=False,
                    **algorithm_options,
                )
                mapped = sample._aggregated_cell_names[factor].map(
                    aggregated_sample.adata.obs[aggregated_key].astype(str)
                )
                sample.adata.obs[key_added] = mapped.where(mapped.notna(), None).astype(object)
                if self._markers_col == key_added:
                    self._markers_col = None
                sample._clustering_factor = factor
                sample.adata.uns["clustering"] = {"algorithm": algorithm, "key": key_added, "factor": factor}
                self._cluster_columns[algorithm] = key_added
                self._active_cluster_col = key_added
                if compute_markers:
                    self.compute_markers()
                with sample.fetch_lock:
                    for cache_key in [k for k in sample._multiscale_cache if key_added in k[2]]:
                        del sample._multiscale_cache[cache_key]
                return key_added

        result_key = self._cluster_on_sample(
            sample,
            algorithm=algorithm,
            resolution=resolution,
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            key_added=key_added,
            compute_markers=compute_markers,
            **algorithm_options,
        )
        self._cluster_columns[algorithm] = result_key
        self._active_cluster_col = result_key
        return result_key

    def compute(self, resolution: float, resolution_level=None, **algorithm_options) -> str:
        """Run clustering (safe in a thread) at the displayed `resolution_level`
        (bins: an aggregation level); marker genes are left for `compute_markers`,
        so the clusters can be shown first. Returns the obs column of the labels.
        """
        options = clustering_kwargs()
        options["resolution"] = resolution
        options.update(algorithm_options)

        return self.run_clustering(
            **options,
            resolution_level=resolution_level,
            compute_markers=False,
        )

    def compute_markers(self, algorithm: str | None = None):
        import scanpy as sc

        cluster_col = self.get_cluster_col(algorithm)
        if cluster_col is None:
            cluster_col = self._active_cluster_col
        if cluster_col is None:
            return

        analysis_sample = self._get_active_analysis_sample()
        adata = analysis_sample.adata_processed
        if adata is None or cluster_col not in adata.obs:
            return

        if self._markers_col == cluster_col:
            return

        report("Ranking marker genes...")

        counts = adata.obs[cluster_col].value_counts()
        groups = [str(label) for label, count in counts.items()
                  if str(label) != "Noise" and count >= 2]
        adata.uns.pop("rank_genes_groups", None)

        if groups and len(counts[counts > 0]) > 1:
            sc.tl.rank_genes_groups(adata, groupby=cluster_col, method=analysis_sample.marker_method,
                                   groups=groups)

        self._markers_col = cluster_col

    def get_marker_genes(self, n_genes=5):
        active_cluster_col = self._active_cluster_col
        analysis_sample = self._get_active_analysis_sample()
        adata_processed = analysis_sample.adata_processed
        if adata_processed is None or 'rank_genes_groups' not in adata_processed.uns or active_cluster_col is None:
            return {}
        if self._markers_col != active_cluster_col:
            return {}

        result = adata_processed.uns['rank_genes_groups']
        clusters = result['names'].dtype.names
        markers = {}

        for cluster in clusters:
            count = min(n_genes, len(result["names"]))
            genes = [result['names'][cluster][i] for i in range(count)]
            scores = [result['scores'][cluster][i] for i in range(count)]
            pvals = [result['pvals_adj'][cluster][i] for i in range(count)]

            # scanpy's log2 fold change of mean normalized expression (cluster vs rest).
            logfcs = [round(float(result['logfoldchanges'][cluster][i]), 2) for i in range(count)]
            
            markers[cluster] = list(zip(
                genes, 
                [round(float(s), 2) for s in scores], 
                [float(p) for p in pvals], 
                [float(l) for l in logfcs]
            ))
            
        return markers

    def _build_cluster_state(self):
        """Derive ordered label list, color map, and active set from obs."""
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None:
            return

        obs = self._base_sample.adata.obs
        if active_cluster_col not in obs.columns:
            return

        # Sort numerically when labels are integers (typical for Leiden)
        raw_labels = obs[active_cluster_col].dropna().astype(str).unique().tolist()
        try:
            self._cluster_labels = sorted(raw_labels, key=lambda x: int(x))
        except ValueError:
            self._cluster_labels = sorted(raw_labels)

        self._active_clusters = set(self._cluster_labels)

        # Assign colors, preserving any colors the user already picked. The map,
        # the DE table and the DE plots all use these (see Component.add_parameter).
        for i, label in enumerate(self._cluster_labels):
            if label not in self._color_map:
                self._color_map[label] = CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]

    def _result_summary(self):
        count = sum(label != "Noise" for label in self._cluster_labels)
        summary = f"{count} clusters found"
        if "Noise" in self._cluster_labels:
            summary += f", {self._count_spots('Noise')} noise observations"
        return summary

    def _count_spots(self, label: str) -> int:
        obs = self._base_sample.adata.obs
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None or active_cluster_col not in obs.columns:
            return 0
        return int((obs[active_cluster_col].astype(str) == label).sum())

    _STEP_DEPENDENTS = {
        "pca": ("neighbors", "tsne", "pca_embedding"),
        "neighbors": ("umap", "leiden"),
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

        report(f"Preprocessing {self.adata.n_obs} observations (normalize, log1p, HVG)...")
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

    # -- views of the state ------------------------------------------------

    def cluster_summary(self, n_genes=settings.MARKER_GENE_COUNT):
        """Clusters of the active result: label, colour, size, shown, marker genes."""
        markers = self.get_marker_genes(n_genes=n_genes)
        return [dict(label=label, colour=self._color_map.get(label, "#cccccc"), size=self._count_spots(label),
                     visible=label in self._active_clusters,
                     markers=[gene for gene, *_ in markers.get(label, [])])
                for label in self._cluster_labels]

    def visible_mask(self, labels, show_not_sequenced=True):
        """Rows whose cluster (in `labels`, aligned to the displayed rows) is shown."""
        labels = pd.Series(labels, dtype=object)
        visible = labels.astype(str).isin(self._active_clusters)
        if show_not_sequenced:
            visible |= labels.isna()
        return np.asarray(visible, dtype=bool)

