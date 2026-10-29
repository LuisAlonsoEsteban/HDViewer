
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
)
from bokeh.layouts import column, row

from .... import settings
from ....settings.config import clustering_kwargs
from ...source.parameter import CategoricalColor
from ...utils.loading import LoadingLogger
from ...utils.palette import CLUSTER_PALETTE, palette_color

NOT_SEQUENCED_VAL = settings.NOT_SEQUENCED_VAL


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

    CLUSTER_COL = settings.CLUSTER_COL

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
        self.cluster_col = self.CLUSTER_COL    # set to the obs column returned by run_clustering

        self._cluster_labels: list[str] = []   # ordered unique cluster ids
        self._color_map: dict[str, str] = {}   # cluster_id → hex color
        self._active_clusters: set[str] = set()

        # Bokeh widget references (populated after clustering)
        self._checkbox_group: CheckboxGroup | None = None
        self._color_pickers: list[ColorPicker] = []
        self._cluster_layout = None             # column() holding cluster widgets

    @property
    def is_clustered(self) -> bool:
        return (
            self.cluster_col in self.component.manager.sample.adata.obs.columns
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

    def compute(self, resolution: float, **algorithm_options) -> str:
        """Run clustering (no Bokeh model updates, safe in a thread); marker genes are
        left for `compute_markers`, so the clusters can be shown first.

        Returns the obs column holding the labels; pass it to `apply_result`.
        """
        manager = self.component.manager

        # Cluster at the level currently displayed (Bins store it as 'leiden_{factor}')
        resolution_level = getattr(manager, "effective_resolution_level", None)

        options = clustering_kwargs()
        options["resolution"] = resolution
        options.update(algorithm_options)

        return manager.sample.run_clustering(
            **options,
            resolution_level=resolution_level,
            compute_markers=False,
        )

    def compute_markers(self):
        """Rank the marker genes of the current clusters (safe in a thread)."""
        self.component.manager.sample.compute_marker_genes()

    def apply_result(self, cluster_col: str, doc, on_done=None):
        """Refresh colors and cluster widgets (must run on the document thread)."""
        self.cluster_col = cluster_col
        self._build_cluster_state()
        self._update_source_colors()
        self._rebuild_cluster_widgets(doc)
        if on_done is not None:
            on_done(doc)

    def get_widgets(self, doc, on_run=None, on_markers=None) -> list:
        """
        Return the full list of widgets for the clustering section:
        resolution slider, run button, status div, reset button, and —
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
            resolution_slider,
            *controls.values(),
            Div(text="Clustering uses expression PCA coordinates. Density methods label outliers as Noise."),
            run_button,
            status.div,
            self._cluster_layout,
        ]

    def _build_cluster_state(self):
        """Derive ordered label list, color map, and active set from obs."""
        obs = self.component.manager.sample.adata.obs
        if self.cluster_col not in obs.columns:
            return

        # Sort numerically when labels are integers (typical for Leiden)
        raw_labels = obs[self.cluster_col].dropna().astype(str).unique().tolist()
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

        if self.cluster_col not in self.component.manager.sample.adata.obs.columns:
            return

        # Re-run de leiden: que los bins agregados no salgan con etiquetas antiguas.
        if hasattr(self.component.manager.sample, 'clear_multiscale_cache'):
            self.component.manager.sample.clear_multiscale_cache()

        self.component.add_parameter(
            attr='color',
            colname=self.cluster_col,
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
            obs_series = self.component.manager.sample.adata.obs[self.cluster_col]
            obs_names = pd.Series(self.component.manager.source.data['obs_name'])
            aligned = obs_names.map(obs_series)
            new_colors = color_param.get_values(aligned.copy())
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            # Columna leiden en el source para el hover ("Cluster: @leiden"),
            # igual que groups.py hace con custom_clusters.
            data[self.cluster_col] = ['—' if pd.isna(c) else str(c) for c in aligned]
            self.component.manager.source.data = data

        self.component.update_hover_tooltips()

    def _apply_visibility_mask(self):
        mgr = self.component.manager
        if self.cluster_col in mgr.source.data and len(mgr.source.data[self.cluster_col]) > 0:
            cluster_ids = pd.Series(mgr.source.data[self.cluster_col])
        elif hasattr(self.component, 'loader'):
            # Bins (lazy): si el color ya no es leiden (p.ej. clusters personalizados),
            # la columna no está en el source. Alinear por obs_name con lo cargado;
            # get_value devolvería el gdf entero y desincroniza view_filter.
            obs_names = pd.Series(mgr.source.data.get('obs_name', []), dtype=object)
            cluster_ids = obs_names.map(mgr.sample.adata.obs[self.cluster_col].astype(str))
        else:
            cluster_ids = mgr.sample.get_value(colname=self.cluster_col)
            
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
        markers = self.component.manager.sample.get_marker_genes(n_genes=settings.MARKER_GENE_COUNT)
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
            cb = CheckboxGroup(labels=[("Noise" if label == "Noise" else "Cluster " + str(label)) + "  (" + str(self._count_spots(label)) + f" {unit})"], active=[0])
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
    # Utilities
    # ------------------------------------------------------------------

    def _result_summary(self):
        count = sum(label != "Noise" for label in self._cluster_labels)
        summary = f"{count} clusters found"
        if "Noise" in self._cluster_labels:
            summary += f", {self._count_spots('Noise')} noise observations"
        return summary

    def _count_spots(self, label: str) -> int:
        obs = self.component.manager.sample.adata.obs
        if self.cluster_col not in obs.columns:
            return 0
        return int((obs[self.cluster_col].astype(str) == label).sum())


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
