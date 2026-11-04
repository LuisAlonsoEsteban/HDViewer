
import pandas as pd
import numpy as np

from bokeh.models import (
    CheckboxGroup,
    ColorPicker,
    Button,
    Div,
    Slider,
    Select,
    Spinner,
    BooleanFilter,
    InlineStyleSheet,
)
from bokeh.layouts import column, row

from .... import settings
from ....settings.config import clustering_kwargs
from ...source.parameter import CategoricalColor
from ...utils.loading import LoadingLogger
from ...utils.palette import CLUSTER_PALETTE, palette_color
from ....startup.utils.progress import report

NOT_SEQUENCED_VAL = settings.NOT_SEQUENCED_VAL

SUPPORTED_ALGORITHMS = (
    "leiden",
    "kmeans",
    "dbscan",
    "hdbscan",
    "gaussian_mixture",
    "agglomerative",
)



def cluster_label_css(color):
    """A cluster's checkbox label in its assigned colour, with a swatch."""
    return (f".bk-input-group label span {{ color: {color}; font-weight: 600; }}"
            f".bk-input-group label span::before {{ content: ''; display: inline-block; width: 10px;"
            f" height: 10px; margin-right: 6px; border-radius: 2px; vertical-align: -1px;"
            f" background: {color}; }}")

class ClusteringManager:
    """
    Manages expression clustering lifecycle for a sample.

    Responsibilities
    ----------------
    - Run/reset clustering via the sample object.
    - Maintain the mapping of cluster label → color.
    - Produce a CheckboxGroup (one label per cluster) that toggles spot
      visibility by masking the ColumnDataSource fed to the glyph renderer.
    - Produce per-cluster ColorPickers that update the color column in the
      same source.

    Usage inside Component
    ----------------------
        self.clustering_manager = ClusteringManager(
            sample=self.sample,
            manager=self.component.manager,   # SourceManager / DynamicSourceManager
        )
        # In get_widgets():
        self.widgets.extend(self.clustering_manager.get_widgets(doc))
    """


    def __init__(self, component):
        """
        Parameters
        ----------
        sample  : object with `.adata`, `.run_clustering(resolution)`,
                  and `.get_colnames()` methods.
        manager : SourceManager / DynamicSourceManager that exposes
                  `.source` (ColumnDataSource) and `.update_colname()`.run
        """
        self.component = component

        self._cluster_columns: dict[str, str] = {}
        self._active_cluster_col: str | None = None
        self._cluster_labels: list[str] = []   # ordered unique cluster ids for the active result
        self._color_map: dict[str, str] = {}  # cluster_id → hex color
        self._active_clusters: set[str] = set()
        self._step_params: dict[str, object] = {}
        self._bound_sample = None

        # Bokeh widget references (populated after clustering)
        self._checkbox_group: CheckboxGroup | None = None
        self._color_pickers: list[ColorPicker] = []
        self._cluster_layout = None             # column() holding cluster widgets

    @property
    def sample(self):
        return self._bound_sample if self._bound_sample is not None else self.component.manager.sample

    @property
    def adata(self):
        return self.sample.adata

    @adata.setter
    def adata(self, value):
        self.sample.adata = value

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
        sample = self.component.manager.sample
        active_col = self._active_cluster_col
        if active_col is None:
            return sample
        factor = self.get_cluster_factor(active_col)
        if factor is not None and hasattr(sample, "_aggregated_samples") and factor in sample._aggregated_samples:
            return sample._aggregated_samples[factor]
        return sample

    # ------------------------------------------------------------------
    # Public lifecycle
    # ------------------------------------------------------------------

    @property
    def is_clustered(self) -> bool:
        return (
            any(column in self.component.manager.sample.adata.obs.columns for column in self._cluster_columns.values())
            and len(self._cluster_labels) > 0
        )

    def run(self, resolution: float, doc, on_done=None):
        """
        Execute clustering (blocking; call inside add_next_tick_callback).

        Parameters
        ----------
        resolution : Leiden resolution passed to sample.run_clustering.
        doc        : Bokeh Document, used for next-tick callbacks.
        on_done    : optional callable(doc) invoked after the source and
                     widgets have been refreshed.
        """
        cluster_col = self.compute(resolution)
        self.compute_markers()
        self.apply_result(cluster_col, doc, on_done)

    def _cluster_on_sample(self, sample, *, algorithm, resolution, n_pcs, n_neighbors,
                          key_added=None, compute_markers=False, **algorithm_options):
        """Run a single clustering algorithm on either a base sample or an aggregated one."""
        previous_sample = self._bound_sample
        previous_step_params = self._step_params.copy()
        self._bound_sample = sample
        self._step_params = {}
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
                # The restored parameters must not claim markers for the old labels.
                if previous_step_params.get("markers") == key_added:
                    previous_step_params.pop("markers")

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
            self._step_params = previous_step_params

    def run_clustering(self, resolution=0.9, n_pcs=30, n_neighbors=15, key_added=None,
                       compute_markers=False, resolution_level=None, *, algorithm="leiden",
                       **algorithm_options):
        """Run clustering and handle bins aggregation in the manager layer."""
        sample = self.component.manager.sample

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

    def compute(self, resolution: float, **algorithm_options) -> str:
        """Run clustering (no Bokeh model updates, safe in a thread); marker genes are
        left for `compute_markers`, so the clusters can be shown first.

        Returns the obs column holding the labels; pass it to `apply_result`.
        """
        manager = self.component.manager

        resolution_level = getattr(manager, "effective_resolution_level", None)

        options = clustering_kwargs()
        options["resolution"] = resolution
        options.update(algorithm_options)

        return self.run_clustering(
            **options,
            resolution_level=resolution_level,
            compute_markers=False,
        )

    # ------------------------------------------------------------------
    # Marker genes
    # ------------------------------------------------------------------

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

        marker_params = self._step_params.get("markers")
        if marker_params == cluster_col:
            return

        report("Ranking marker genes...")

        counts = adata.obs[cluster_col].value_counts()
        groups = [str(label) for label, count in counts.items()
                  if str(label) != "Noise" and count >= 2]
        adata.uns.pop("rank_genes_groups", None)

        if groups and len(counts[counts > 0]) > 1:
            sc.tl.rank_genes_groups(adata, groupby=cluster_col, method=analysis_sample.marker_method,
                                   groups=groups)

        self._step_params["markers"] = cluster_col

    def get_marker_genes(self, n_genes=5):
        active_cluster_col = self._active_cluster_col
        analysis_sample = self._get_active_analysis_sample()
        adata_processed = analysis_sample.adata_processed
        if adata_processed is None or 'rank_genes_groups' not in adata_processed.uns or active_cluster_col is None:
            return {}
        if self._step_params.get("markers") != active_cluster_col:
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

    # ------------------------------------------------------------------
    # Result application and widgets
    # ------------------------------------------------------------------

    def apply_result(self, cluster_col: str, doc, on_done=None):
        """Refresh colors and cluster widgets (must run on the document thread)."""
        self._active_cluster_col = cluster_col
        self._build_cluster_state()
        self._update_source_colors()
        self._rebuild_cluster_widgets(doc)
        if on_done is not None:
            on_done(doc)

    def get_widgets(self, doc, on_run=None, on_markers=None) -> list:
        """
        Return the full list of widgets for the clustering section:
        algorithm selector, resolution slider, run button, status div, and —
        after clustering — the cluster checkbox/color widgets.

        Safe to call multiple times; widget state is preserved.

        on_run(doc) is called once the clusters are shown, on_markers(doc) once their
        marker genes are ranked (that runs afterwards, in the background).
        """
        resolution_slider = Slider(
            start=settings.CLUSTER_SLIDER_MIN, end=settings.CLUSTER_SLIDER_MAX,
            value=settings.CLUSTER_SLIDER_DEFAULT, step=settings.CLUSTER_SLIDER_STEP,
            title="Clustering resolution",
        )
        if not hasattr(self, "_algorithm_controls"):
            defaults = clustering_kwargs()

            self._algorithm_controls = {
                "algorithm": Select(title="Clustering algorithm", value=defaults["algorithm"], options=[
                    ("leiden", "Leiden"), ("kmeans", "K-means"), ("dbscan", "DBSCAN"),
                    ("hdbscan", "HDBSCAN"), ("gaussian_mixture", "Gaussian mixtures"),
                    ("agglomerative", "Agglomerative")]),
                "n_clusters": Spinner(title="Number of clusters", low=1, step=1, value=defaults["n_clusters"]),
                "eps": Spinner(title="DBSCAN epsilon (PCA distance)", low=0.000001, step=0.1, value=defaults["eps"]),
                "min_samples": Spinner(title="Minimum samples", low=1, step=1, value=defaults["min_samples"]),
                "min_cluster_size": Spinner(title="Minimum cluster size", low=2, step=1, value=defaults["min_cluster_size"]),
                "covariance_type": Select(title="Covariance type", value=defaults["covariance_type"],
                    options=["full", "tied", "diag", "spherical"]),
                "linkage": Select(title="Linkage", value=defaults["linkage"],
                    options=["ward", "complete", "average", "single"]),
            }

        controls = self._algorithm_controls

        # A rebuilt panel gets its own resolution widget; replace the old handler.
        algorithm_select = controls["algorithm"]
        if hasattr(self, "_algorithm_visibility_callback"):
            algorithm_select.remove_on_change("value", self._algorithm_visibility_callback)

        def _show_algorithm_options(attr, old, new):
            algorithm = algorithm_select.value
            resolution_slider.visible = algorithm == "leiden"
            fields = {
                "n_clusters": ("kmeans", "gaussian_mixture", "agglomerative"),
                "eps": ("dbscan",), "min_samples": ("dbscan", "hdbscan"),
                "min_cluster_size": ("hdbscan",), "covariance_type": ("gaussian_mixture",),
                "linkage": ("agglomerative",),
            }
            for key, algorithms in fields.items():
                controls[key].visible = algorithm in algorithms

        self._algorithm_visibility_callback = _show_algorithm_options
        algorithm_select.on_change("value", _show_algorithm_options)
        _show_algorithm_options(None, None, None)

        run_button = Button(label="Run clustering", button_type="success")

        status = LoadingLogger(doc)
        status.set_message("Run clustering to see results")

        # Lazy-build cluster widget placeholder so we can update it in-place
        if self._cluster_layout is None:
            self._cluster_layout = column(
                name="cluster_controls",
                children=self._build_cluster_control_widgets(doc),
            )

        def _on_run():
            if status.running:
                return

            resolution = resolution_slider.value
            algorithm_options = {key: widget.value for key, widget in controls.items()}
            integer_fields = ("n_clusters", "min_samples", "min_cluster_size")

            if any(algorithm_options[key] is None or algorithm_options[key] != int(algorithm_options[key])
                   for key in integer_fields):
                status.error("Cluster counts and minimum samples must be whole numbers.")
                return

            for key in integer_fields:
                algorithm_options[key] = int(algorithm_options[key])

            def _on_markers_done(_):
                run_button.disabled = False
                self._markers_pending = False
                self._rebuild_cluster_widgets(doc)
                status.success(f"{self._result_summary()}, marker genes ready")
                if on_markers is not None:
                    on_markers(doc)

            def _on_markers_error(exc):
                run_button.disabled = False
                self._markers_pending = False
                self._rebuild_cluster_widgets(doc)
                status.error(f"Marker genes failed: {str(exc)[:60]}")

            def _on_done(cluster_col):
                try:
                    # Show the clusters now; the marker genes follow in the background
                    self._markers_pending = True
                    self.apply_result(cluster_col, doc, on_done=on_run)
                except Exception as exc:
                    run_button.disabled = False
                    self._markers_pending = False
                    status.error(str(exc)[:60])
                    return

                status.run(f"{self._result_summary()} · marker genes",
                           self.compute_markers,
                           on_done=_on_markers_done, on_error=_on_markers_error)

            def _on_error(exc):
                run_button.disabled = False
                status.error(str(exc)[:60])

            # Stays disabled until the marker genes are ready too
            run_button.disabled = True
            status.run(
                f"Clustering ({algorithm_options['algorithm']})",
                lambda: self.compute(resolution, **algorithm_options),
                on_done=_on_done,
                on_error=_on_error,
            )

        run_button.on_click(_on_run)

        return [
            algorithm_select,
            resolution_slider,
            *(widget for key, widget in controls.items() if key != "algorithm"),
            Div(text="Clustering uses expression PCA coordinates. Density methods label outliers as Noise."),
            run_button,
            status.div,
            self._cluster_layout,
        ]

    def _build_cluster_state(self):
        """Derive ordered label list, color map, and active set from obs."""
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None:
            return

        obs = self.component.manager.sample.adata.obs
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
                self._color_map[label] = palette_color(CLUSTER_PALETTE, i)

    def _update_source_colors(self):
        """
        Write a 'cluster_color' column into the ColumnDataSource so the
        glyph renderer can use field('cluster_color') for fill_color.
        """
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None or active_cluster_col not in self.component.manager.sample.adata.obs.columns:
            return

        # Re-run de leiden: que los bins agregados no salgan con etiquetas antiguas.
        if hasattr(self.component.manager.sample, 'clear_multiscale_cache'):
            self.component.manager.sample.clear_multiscale_cache()

        self.component.add_parameter(
            attr='color',
            colname=active_cluster_col,
            static=False,
            not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.component.alias.lower()}s"]["color"],
            colname_hovered=True,
        )

        self.component.manager.update_param('alpha', avoid_scaling=True)

        # Bins (lazy loader): recargar por tiles. NO inyectar el source a mano,
        # eso desincroniza source.data con self.indexes y revienta remove_further_data.
        if hasattr(self.component, 'loader'):
            self.component.loader.force_reload()
            self.component.update_hover_tooltips()
            return


        # Cells: la inyección manual no cuela el cluster (obs_name vacío / no casa).
        # Re-fetch del viewport -> pasa por get_value y alinea adata.obs con el gdf,
        # igual que hace apply_colors_to_column para componentes sin loader.
        if self.component.alias.lower() == 'cell':
            # Cells son TILED (manager.source vacío). Igual que groups:
            # el color param ya está puesto y leiden está en adata.obs;
            # solo forzamos el re-render de los tiles vía change_glyph_view.
            gdf = self.component.manager.sample.gdf
            self.component._view = [True] * len(gdf)
            self.component.change_glyph_view()
            self.component.update_hover_tooltips()
            return

        color_param = self.component.manager.get_param('color')
        if color_param is not None:
            obs_series = self.component.manager.sample.adata.obs[active_cluster_col]
            obs_names = pd.Series(self.component.manager.source.data['obs_name'])
            aligned = obs_names.map(obs_series)
            new_colors = color_param.get_values(aligned.copy())
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            # Columna del algoritmo actual en el source para el hover.
            data[active_cluster_col] = ['—' if pd.isna(c) else str(c) for c in aligned]
            self.component.manager.source.data = data

        self.component.update_hover_tooltips()

    def _apply_visibility_mask(self):
        mgr = self.component.manager
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None:
            return
        if active_cluster_col in mgr.source.data and len(mgr.source.data[active_cluster_col]) > 0:
            cluster_ids = pd.Series(mgr.source.data[active_cluster_col])
        elif hasattr(self.component, 'loader'):
            # Bins (lazy): si el color ya no es leiden (p.ej. clusters personalizados),
            # la columna no está en el source. Alinear por obs_name con lo cargado;
            # get_value devolvería el gdf entero y desincroniza view_filter.
            obs_names = pd.Series(mgr.source.data.get('obs_name', []), dtype=object)
            cluster_ids = obs_names.map(mgr.sample.adata.obs[active_cluster_col].astype(str))
        else:
            cluster_ids = mgr.sample.get_value(colname=active_cluster_col)
            
        visible = cluster_ids.astype(str).isin(self._active_clusters)

        if self.component.show_not_sequenced:
            visible = visible | cluster_ids.isna()
        visible = np.asarray(visible.tolist(), dtype=bool)
        
        if hasattr(mgr, "view_filter"):
            mgr.view_filter = visible.copy()

        self.component._view = visible.tolist()
        self.component.change_glyph_view()

    # ------------------------------------------------------------------
    # Widget builders
    # ------------------------------------------------------------------

    def _build_cluster_control_widgets(self, doc) -> list:
        if not self.is_clustered:
            return []
            
        from bokeh.models import CheckboxGroup, Button
        from bokeh.models import Div as _Div

        markers = self.get_marker_genes(n_genes=settings.MARKER_GENE_COUNT)
        self._cluster_checkboxes = []

        widgets = []
        select_all_btn = Button(label="Select all / Deselect all", button_type="default", name="cluster_select_all")
        widgets.append(select_all_btn)

        def _on_any_change(attr, old, new):
            if getattr(self, "_bulk_update", False):
                return
            active = set()
            for cb, lbl in zip(self._cluster_checkboxes, self._cluster_labels):
                if 0 in cb.active:
                    active.add(lbl)
            self._active_clusters = active
            doc.add_next_tick_callback(self._apply_visibility_mask)

        for label in self._cluster_labels:
            unit = f"{self.component.alias.lower()}s"  # spots / bins / cells
            color = self._color_map.get(str(label), "#cccccc")
            cb = CheckboxGroup(labels=[("Noise" if label == "Noise" else "Cluster " + str(label)) + "  (" + str(self._count_spots(label)) + f" {unit})"], active=[0],
                               stylesheets=[InlineStyleSheet(css=cluster_label_css(color))])
            cb.on_change("active", _on_any_change)
            self._cluster_checkboxes.append(cb)
            entry = markers.get(str(label)) or markers.get(label)
            if entry:
                gns = ", ".join(str(t[0]) for t in entry)
            elif getattr(self, "_markers_pending", False):
                gns = "computing marker genes…"
            else:
                gns = "-"
            gdiv = _Div(text="&nbsp;&nbsp;&nbsp;&nbsp;<span style='color:#aaaaaa;font-size:11px'>" + gns + "</span>", width=270)
            widgets.append(cb)
            widgets.append(gdiv)

        self._active_clusters = set(self._cluster_labels)
        self._apply_visibility_mask()

        def _on_select_all():
            all_on = all(0 in cb.active for cb in self._cluster_checkboxes)
            self._bulk_update = True
            for cb in self._cluster_checkboxes:
                cb.active = [] if all_on else [0]
            self._bulk_update = False
            self._active_clusters = set() if all_on else set(self._cluster_labels)
            doc.add_next_tick_callback(self._apply_visibility_mask)
        select_all_btn.on_click(_on_select_all)

        return widgets

    def _rebuild_cluster_widgets(self, doc):
        """Replace the children of _cluster_layout in-place."""
        if self._cluster_layout is None:
            return

        def _update():
            self._cluster_layout.children = (
                self._build_cluster_control_widgets(doc)
            )

        doc.add_next_tick_callback(_update)

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------

    def _result_summary(self):
        count = sum(label != "Noise" for label in self._cluster_labels)
        summary = f"{count} clusters found"
        if "Noise" in self._cluster_labels:
            summary += f", {self._count_spots('Noise')} noise observations"
        return summary

    def _count_spots(self, label: str) -> int:
        obs = self.component.manager.sample.adata.obs
        active_cluster_col = self._active_cluster_col
        if active_cluster_col is None or active_cluster_col not in obs.columns:
            return 0
        return int((obs[active_cluster_col].astype(str) == label).sum())


    def apply_colors_to_column(self, colname: str):
        self.component.add_parameter(
            attr='color',
            colname=colname,
            static=False,
            not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.component.alias.lower()}s"]["color"],
            colname_hovered=True,
        )
        self.component.manager.update_param('alpha', avoid_scaling=True)
        if not hasattr(self.component, 'loader'):
            self.component.manager.fetch_data()
        self.component.update_hover_tooltips()

    # ------------------------------------------------------------------
    # Preprocessing and caching
    # ------------------------------------------------------------------

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
