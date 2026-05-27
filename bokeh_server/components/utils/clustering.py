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
                  `.source` (ColumnDataSource) and `.update_colname()`.
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

    def reset(self, doc, on_done=None):
        """Remove clustering results and tear down cluster widgets."""
        obs = self.component.manager.sample.adata.obs
        if self.CLUSTER_COL in obs.columns:
            obs.drop(columns=[self.CLUSTER_COL], inplace=True)

        self._cluster_labels = []
        self._color_map = {}
        self._active_clusters = set()
        self._checkbox_group = None
        self._color_pickers = []

        self._clear_source_colors()
        self._rebuild_cluster_widgets(doc)   # renders empty placeholder

        if on_done is not None:
            on_done(doc)

    def get_widgets(self, doc) -> list:
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
        reset_button = Button(label="Reset clustering", button_type="warning")
        reset_button.visible = self.is_clustered

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
                    )
                    status_div.text = (
                        f"<p style='color:green;font-size:12px'>"
                        f"{len(self._cluster_labels)} clusters found</p>"
                    )
                    reset_button.visible = True
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

        def _on_reset():
            def _do():
                self.reset(doc=doc)
                status_div.text = (
                    "<p style='color:gray;font-size:12px'>"
                    "Run clustering to see results</p>"
                )
                reset_button.visible = False

            doc.add_next_tick_callback(_do)

        run_button.on_click(_on_run)
        reset_button.on_click(_on_reset)

        return [
            resolution_slider,
            run_button,
            status_div,
            reset_button,
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
            not_sequenced_val="#000",
            colname_hovered=True,
        )

        self.component.manager.update_param('alpha', avoid_scaling=True)

    def _clear_source_colors(self):
        """Remove cluster columns from the source after a reset."""
        # self.component.
        pass

    def _apply_visibility_mask(self):
        cluster_ids = self.component.manager.sample.get_value(colname=self.CLUSTER_COL).values
        visible = [
            (
                c in self._active_clusters 
                or
                c not in self._cluster_labels # NaN values
            )
            for c in cluster_ids
        ]
        self.component.glyph_renderer.view.filters = [BooleanFilter(booleans=visible)]

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

        widgets = []

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

        def _on_checkbox_change(attr, old, new):
            self._active_clusters = {
                self._cluster_labels[i] for i in new
            }

            def _patch():
                self._apply_visibility_mask()

            doc.add_next_tick_callback(_patch)

        self._checkbox_group.on_change("active", _on_checkbox_change)
        widgets.append(self._checkbox_group)

        # --- ColorPicker: one per cluster ------------------------------
        self._color_pickers = []
        color_picker_rows = []

        for label in self._cluster_labels:
            picker = ColorPicker(
                title=f"Cluster {label}",
                color=self._color_map[label],
                width=120,
            )

            def _on_color_change(attr, old, new, _label=label):
                self._color_map[_label] = new

                def _repaint():
                    self._update_source_colors()
                    self._apply_visibility_mask()

                doc.add_next_tick_callback(_repaint)

            picker.on_change("color", _on_color_change)
            self._color_pickers.append(picker)

        # Group pickers in rows of 3 to avoid an endless vertical list
        for i in range(0, len(self._color_pickers), 3):
            color_picker_rows.append(
                row(*self._color_pickers[i:i + 3])
            )

        widgets.extend(color_picker_rows)
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
        if self.CLUSTER_COL not in obs.columns:
            return 0
        return int((obs[self.CLUSTER_COL].astype(str) == label).sum())