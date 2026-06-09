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
            self.CLUSTER_COL in self.component.manager.sample.adata.obs.columns
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
        self.component.manager.sample.run_clustering(resolution=resolution)
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
        if self.CLUSTER_COL not in obs.columns:
            return

        # Sort numerically when labels are integers (typical for Leiden)
        raw_labels = obs[self.CLUSTER_COL].astype(str).unique().tolist()
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

        if self.CLUSTER_COL not in self.component.manager.sample.adata.obs.columns:
            return

        self.component.add_parameter(
            attr='color',
            colname=self.CLUSTER_COL,
            static=False,
            not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.component.alias.lower()}s"]["color"],
            colname_hovered=False,
        )

        self.component.manager.update_param('alpha', avoid_scaling=True)

    def _apply_visibility_mask(self):
        cluster_ids = self.component.manager.sample.get_value(
            colname=self.CLUSTER_COL
        )

        visible = cluster_ids.isin(self._active_clusters)

        if 0 in self.component.notsequenced_checkbox.active:
            visible |= cluster_ids.isna()

        self.component._view = visible.to_list()
        self.component.change_glyph_view()

    # ------------------------------------------------------------------
    # Widget builders
    # ------------------------------------------------------------------

    def _build_cluster_control_widgets(self, doc) -> list:
        """
        Return a list of widgets for the current cluster state.
        Empty list when no clustering has been run yet.
        """
        if not self.is_clustered:
            return []

        # --- CheckboxGroup: one entry per cluster ----------------------
        labels = [
            f"Cluster {label}  ({self._count_spots(label)} spots)"
            for label in self._cluster_labels
        ]
        self._checkbox_group = CheckboxGroup(
            labels=labels,
            active=list(range(len(self._cluster_labels))),
            name="cluster_checkbox_group",
        )

        self._apply_visibility_mask() # Make all spots appear

        def _on_checkbox_change(attr, old, new):
            self._active_clusters = {
                self._cluster_labels[i] for i in new
            }

            def _patch():
                self._apply_visibility_mask()

            doc.add_next_tick_callback(_patch)

        self._checkbox_group.on_change("active", _on_checkbox_change)

        return [self._checkbox_group]

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
        if self.CLUSTER_COL not in obs.columns:
            return 0
        return int((obs[self.CLUSTER_COL].astype(str) == label).sum())


    #def apply_colors_to_column(self, colname: str):
    #    """
    #    Aplica colores a cualquier columna categórica de adata.obs.
    #    Mismo flujo que _update_source_colors pero con colname arbitrario.
    #    """
    #    from .. import NOT_SEQUENCED_VAL
#
    #    self.component.add_parameter(
    #        attr='color',
    #        colname=colname,
    #        static=False,
    #        not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.component.alias.lower()}s"]["color"],
    #        colname_hovered=True,
    #    )
    #    self.component.manager.update_param('alpha', avoid_scaling=True)

    def apply_colors_to_column(self, colname: str):
        from .. import NOT_SEQUENCED_VAL
        self.component.add_parameter(
            attr='color',
            colname=colname,
            static=False,
            not_sequenced_val=NOT_SEQUENCED_VAL[f"{self.component.alias.lower()}s"]["color"],
            colname_hovered=True,
        )
        self.component.manager.update_param('alpha', avoid_scaling=True)