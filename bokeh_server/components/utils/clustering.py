import pandas as pd
import numpy as np
from bokeh.models import (
    CheckboxGroup,
    ColorPicker,
    Button,
    Div,
    Slider,

    BooleanFilter,
)
from bokeh.layouts import column, row
from bokeh.palettes import Category20, Category10

from .. import NOT_SEQUENCED_VAL
from ...parameter import CategoricalColor

CLUSTER_PALETTES = [*Category20[20], *Category10[10]]

class ClusteringManager:
    """
    Manages Leiden clustering lifecycle for a sample.

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

    CLUSTER_COL = "leiden"

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
        manager = self.component.manager
        # Cluster at the level currently displayed (Bins store it as 'leiden_{factor}')
        resolution_level = getattr(manager, "effective_resolution_level", None)
        self.cluster_col = manager.sample.run_clustering(
            resolution=resolution,
            resolution_level=resolution_level,
        )
        self._build_cluster_state()
        self._update_source_colors()
        self._rebuild_cluster_widgets(doc)
        if on_done is not None:
            on_done(doc)

    def get_widgets(self, doc, on_run=None) -> list:
        """
        Return the full list of widgets for the clustering section:
        resolution slider, run button, status div, reset button, and —
        after clustering — the cluster checkbox/color widgets.

        Safe to call multiple times; widget state is preserved.
        """
        resolution_slider = Slider(
            start=0.1, end=2.0, value=1.0, step=0.1,
            title="Clustering resolution",
        )
        run_button = Button(label="Run clustering", button_type="success")
        status_div = Div(
            text="<p style='color:gray;font-size:12px'>"
                 "Run clustering to see results</p>"
        )

        # Lazy-build cluster widget placeholder so we can update it in-place
        if self._cluster_layout is None:
            self._cluster_layout = column(
                name="cluster_controls",
                children=self._build_cluster_control_widgets(doc),
            )

        def _on_run():
            def _do():
                run_button.disabled = True
                status_div.text = (
                    "<p style='color:gray;font-size:12px'>Clustering…</p>"
                )
                try:
                    self.run(
                        resolution=resolution_slider.value,
                        doc=doc,
                        on_done=on_run
                    )
                    status_div.text = (
                        f"<p style='color:green;font-size:12px'>"
                        f"{len(self._cluster_labels)} clusters found</p>"
                    )
                except Exception as exc:
                    import traceback
                    print(traceback.format_exc())
                    status_div.text = (
                        f"<p style='color:red;font-size:12px'>"
                        f"Error: {str(exc)[:60]}</p>"
                    )
                finally:
                    run_button.disabled = False

            doc.add_next_tick_callback(_do)

        run_button.on_click(_on_run)

        return [
            resolution_slider,
            run_button,
            status_div,
            self._cluster_layout,
        ]

    def _build_cluster_state(self):
        """Derive ordered label list, color map, and active set from obs."""
        obs = self.component.manager.sample.adata.obs
        if self.cluster_col not in obs.columns:
            return

        # Sort numerically when labels are integers (typical for Leiden)
        raw_labels = obs[self.cluster_col].astype(str).unique().tolist()
        try:
            self._cluster_labels = sorted(raw_labels, key=lambda x: int(x))
        except ValueError:
            self._cluster_labels = sorted(raw_labels)

        self._active_clusters = set(self._cluster_labels)

        # Assign colors, preserving any colors the user already picked
        palette = CLUSTER_PALETTES
        for i, label in enumerate(self._cluster_labels):
            if label not in self._color_map:
                self._color_map[label] = palette[i % len(palette)]

    def _update_source_colors(self):
        """
        Write a 'cluster_color' column into the ColumnDataSource so the
        glyph renderer can use field('cluster_color') for fill_color.
        """

        if self.cluster_col not in self.component.manager.sample.adata.obs.columns:
            return

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

        print(f"[CLUSTER COLORS] alias={self.component.alias} "
              f"loader={hasattr(self.component, 'loader')} "
              f"obs_in_source={len(self.component.manager.source.data.get('obs_name', []))}")

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
            obs_series = self.component.manager.sample.adata.obs[self.CLUSTER_COL]
            obs_names = pd.Series(self.component.manager.source.data['obs_name'])
            aligned = obs_names.map(obs_series)
            new_colors = color_param.get_values(aligned.copy())
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            self.component.manager.source.data = data

        self.component.update_hover_tooltips()

    def _apply_visibility_mask(self):
        mgr = self.component.manager
        if self.CLUSTER_COL in mgr.source.data and len(mgr.source.data[self.CLUSTER_COL]) > 0:
            cluster_ids = pd.Series(mgr.source.data[self.CLUSTER_COL])
        else:
            cluster_ids = mgr.sample.get_value(colname=self.CLUSTER_COL)
            
        cluster_num = pd.to_numeric(cluster_ids, errors="coerce")
        active_num = set(pd.to_numeric(pd.Series(list(self._active_clusters)), errors="coerce").dropna())
        visible = cluster_num.isin(active_num)

        if 0 in self.component.notsequenced_checkbox.active:
            visible = visible | cluster_num.isna()
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
        markers = self.component.manager.sample.get_marker_genes(n_genes=5)
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
            cb = CheckboxGroup(labels=["Cluster " + str(label) + "  (" + str(self._count_spots(label)) + " spots)"], active=[0])
            cb.on_change("active", _on_any_change)
            self._cluster_checkboxes.append(cb)
            entry = markers.get(str(label)) or markers.get(label)
            gns = ", ".join(str(t[0]) for t in entry) if entry else "-"
            gdiv = _Div(text="&nbsp;&nbsp;&nbsp;&nbsp;<span style='color:#aaaaaa;font-size:11px'>" + gns + "</span>", width=380)
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
