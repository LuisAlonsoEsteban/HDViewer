
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
from ...source.parameter import CategoricalColor
from ...utils.loading import LoadingLogger
from ....settings.config import clustering_kwargs
from ....core.clustering import ClusteringState, SUPPORTED_ALGORITHMS  # noqa: F401  (re-exported)

NOT_SEQUENCED_VAL = settings.NOT_SEQUENCED_VAL




def cluster_label_css(color):
    """A cluster's checkbox label in its assigned colour, with a swatch."""
    return (f".bk-input-group label span {{ color: {color}; font-weight: 600; }}"
            f".bk-input-group label span::before {{ content: ''; display: inline-block; width: 10px;"
            f" height: 10px; margin-right: 6px; border-radius: 2px; vertical-align: -1px;"
            f" background: {color}; }}")

class ClusteringManager(ClusteringState):
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
        super().__init__()

        # Bokeh widget references (populated after clustering)
        self._checkbox_group: CheckboxGroup | None = None
        self._color_pickers: list[ColorPicker] = []
        self._cluster_layout = None             # column() holding cluster widgets

















    # ------------------------------------------------------------------
    # Public lifecycle
    # ------------------------------------------------------------------

    @property
    def _base_sample(self):
        return self.component.manager.sample

    def compute(self, resolution: float, **algorithm_options) -> str:
        """Cluster at the layer's displayed level (see ClusteringState.compute)."""
        level = getattr(self.component.manager, "effective_resolution_level", None)
        return super().compute(resolution, resolution_level=level, **algorithm_options)


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




    # ------------------------------------------------------------------
    # Marker genes
    # ------------------------------------------------------------------



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
            self.component._view = np.ones(len(gdf), dtype=bool)
            self.component.change_glyph_view()
            self.component.update_hover_tooltips()
            return

        color_param = self.component.manager.get_param('color')
        if color_param is not None:
            obs_series = self.component.manager.sample.adata.obs[active_cluster_col]
            obs_names = pd.Series(self.component.manager.loaded_obs_names)
            aligned = obs_names.map(obs_series)
            new_colors = color_param.get_values(aligned.copy())
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            # Columna del algoritmo actual en el source para el hover.
            data[active_cluster_col] = aligned.astype(str).where(aligned.notna(), '—').tolist()
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
            obs_names = pd.Series(mgr.loaded_obs_names, dtype=object)
            cluster_ids = obs_names.map(mgr.sample.adata.obs[active_cluster_col].astype(str))
        else:
            cluster_ids = mgr.sample.get_value(colname=active_cluster_col)
            
        visible = cluster_ids.astype(str).isin(self._active_clusters)

        if self.component.show_not_sequenced:
            visible = visible | cluster_ids.isna()
        visible = np.asarray(visible, dtype=bool)

        if hasattr(mgr, "view_filter"):
            mgr.view_filter = visible.copy()

        self.component._view = visible
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







