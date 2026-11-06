
import numpy as np
import pandas as pd

from bokeh.models import Button, CheckboxGroup, Div, RadioButtonGroup, Select
from bokeh.layouts import row

from .... import settings
from ...utils.loading import LoadingLogger

from ....core.de_sets import cluster_label, de_colours, de_masks, de_options
from ...operations.layer.de import (
    DEAnalysis,
    build_de_plots,
    DEPresenter
)


class DEManager:
    """Container for DE UI widgets and the analysis lifecycle.

    The class intentionally mirrors the current behavior that was embedded in the
    groups manager so it can be moved without changing the surrounding app
    structure yet.
    """

    DE_SOURCES = ("groups", "clusters")
    CLUSTER_PREFIX = settings.DE_CLUSTER_PREFIX

    def __init__(self, component):
        self.component = component
        self.sample = self._resolve_sample()
        self.de_analysis = DEAnalysis(self.sample)

        if hasattr(self.component, "groups_manager"):
            self.component.groups_manager.add_change_listener(self._groups_changed)

        self._color_map = {}
        self._de_groups_hint = None
        self._de_groups_checkbox = None
        self._de_clusters_hint = None
        self._de_clusters_checkbox = None
        self._de_selected = set()
        self._de_selected_clusters = None
        self._de_results = {}
        self._de_merged = None
        self._de_clusters = None
        self._de_results_source = None
        self._de_generation = 0
        self._de_status_div = None
        self._de_logger = None

    def _resolve_sample(self):
        manager = getattr(self.component, "manager", None)
        if manager is not None and hasattr(manager, "sample"):
            return manager.sample
        return getattr(self.component, "sample", None)

    @property
    def groups(self):
        if hasattr(self.component, "groups_manager"):
            return self.component.groups_manager
        return None

    def _de_options(self, source):
        """Names of the sets available for DE, in display order."""
        return de_options(source, getattr(self.component, "groups_manager", None), self._clustering())

    def _cluster_label(self, name):
        return cluster_label(name)

    def _clustering(self):
        manager = getattr(self.component, "clustering_manager", None)
        return manager if manager is not None and getattr(manager, "is_clustered", False) else None

    def _de_chosen(self, source):
        if source == "groups":
            return self._de_selected
        return self._de_selected_clusters or set()

    def _de_masks(self, source, names):
        """Boolean membership per set, aligned with the sample's observations."""
        return de_masks(self.sample, source, names, getattr(self.component, "groups_manager", None),
                        self._clustering())

    def _de_colors(self, source):
        return de_colours(source, getattr(self.component, "groups_manager", None), self._clustering())

    def _get_de_logger(self, doc):
        if getattr(self, "_de_logger", None) is None:
            self._de_logger = LoadingLogger(doc, width=270, font_size=11)
            self._de_status_div = self._de_logger.div
        return self._de_logger

    def _displayed_factor(self):
        manager = getattr(self.component, "manager", None)
        if manager is None or self.sample is None:
            return 1
        return self.sample.get_aggregation_factor(getattr(manager, "effective_resolution_level", None))

    def _de_level_select(self):
        factors = getattr(self.sample, "AGGREGATION_FACTORS", None)
        if not factors:
            return None

        def displayed_label():
            factor = self._displayed_factor()
            return f"Displayed level ({'native' if factor == 1 else f'{factor} × {factor}'})"

        select = Select(
            title="Resolution",
            value="displayed",
            options=[("displayed", displayed_label()), ("1", "Native bins")]
            + [(str(f), f"{f} × {f} bins (summed counts)") for f in factors if f > 1],
            width=235,
        )

        def _on_context_change():
            select.options = [("displayed", displayed_label())] + list(select.options[1:])

        on_context_change = getattr(self.component.manager, "on_context_change", None)
        if on_context_change is not None:
            on_context_change(_on_context_change)

        return select

    def _de_factor_for(self, level_select):
        if level_select is None:
            return 1
        return self._displayed_factor() if level_select.value == "displayed" else int(level_select.value)

    def _run_de_thread(self, active, mode, doc, masks=None, colors=None, source="groups", factor=1):
        """Run the independent analysis and display its result."""
        logger = self._get_de_logger(doc)
        if logger.running:
            return

        masks = self._de_masks("groups", active) if masks is None else masks
        masks = {name: masks[name] for name in active if name in masks}
        colors = dict(self._color_map) if colors is None else dict(colors)
        self._de_clusters = list(active)
        self._de_results_source = source
        generation = self._de_generation
        unit, units = {"groups": ("group", "groups"), "clusters": ("cluster", "clusters")}.get(
            source, ("set", "groups/clusters")
        )

        def work():
            result = self.de_analysis.run(masks, mode, source=source, factor=factor)
            presenter = DEPresenter(
                alias=self.component.alias,
                sample_type=getattr(self.component, "sample_type", None),
                colors=colors,
            )
            return result, presenter.prepare(result)

        def on_done(output):
            # Report presentation failures too: otherwise the spinner never clears.
            try:
                show(output)
            except Exception as exc:
                logger.error(exc)

        def show(output):
            if generation != self._de_generation:
                logger.set_message("DE results discarded because the comparison groups changed.")
                return
            result, view = output
            self._de_header = view.header
            self._de_results = result.results
            self._de_merged = view.merged
            de_layout = getattr(self.component, "_figure_de_layout", None)
            if de_layout is not None:
                plots = build_de_plots(view.plot_data)
                de_layout._comparison_owner = self
                de_layout.children = [Div(text=view.table_html, sizing_mode="stretch_width"), plots]
            obs_word = f"{self.component.alias.lower()}s" if factor <= 1 else f"{factor} × {factor} {self.component.alias.lower()}s"
            left_out = f", {result.n_overlap} overlapping {obs_word} left out" if result.n_overlap else ""
            count = len(result.results)
            logger.success(f"DE done: {count} {unit if count == 1 else units}{left_out}, results below the map")

        hint = f"{units} vs rest" if mode == 0 else f"between selected {units}"
        logger.run(f"DE ({hint})", work, on_done=on_done, on_error=logger.error)

    def _groups_changed(self):
        self._refresh_de_groups()
        self._clear_de_if_cluster_deleted()

    def _refresh_de_groups(self):
        """Keep both DE lists in step with the groups and clusters that exist."""
        if self._de_groups_checkbox is None:
            return

        groups = self._de_options("groups")
        self._de_selected &= set(groups)
        self._de_groups_checkbox.labels = groups
        self._de_groups_checkbox.active = [i for i, n in enumerate(groups) if n in self._de_selected]
        self._de_groups_hint.text = (
            "" if groups else "<p style='color:#888;font-size:11px'>Create groups in the Groups tab.</p>"
        )

        clusters = self._de_options("clusters")
        self._de_selected_clusters = (self._de_selected_clusters or set()) & set(clusters)
        self._de_clusters_checkbox.labels = clusters
        self._de_clusters_checkbox.active = [
            i for i, n in enumerate(clusters) if n in self._de_selected_clusters
        ]
        self._de_clusters_hint.text = (
            "" if clusters else "<p style='color:#888;font-size:11px'>Run clustering in the Clustering tab.</p>"
        )

    def refresh_de_clusters(self):
        """Called after (re)clustering: the cluster list and any cluster DE are stale."""
        self._de_selected_clusters = set()
        if getattr(self, "_de_results_source", None) in ("clusters", "mixed"):
            self._clear_de_results()
        self._refresh_de_groups()

    def _clear_de_results(self):
        self._de_generation += 1
        self._de_clusters = None
        self._de_results = {}
        de_layout = getattr(self.component, "_figure_de_layout", None)
        if de_layout is not None and getattr(de_layout, '_comparison_owner', self) in (None, self):
            de_layout.children = []
            de_layout._comparison_owner = None
        if self._de_status_div is not None:
            self._de_status_div.text = ""

    def get_widgets(self, doc):
        """Build the DE controls panel."""
        self._de_groups_hint = Div(width=270)
        self._de_groups_checkbox = CheckboxGroup(labels=[], active=[])
        self._de_clusters_hint = Div(width=270)
        self._de_clusters_checkbox = CheckboxGroup(labels=[], active=[])
        self._refresh_de_groups()

        def _on_groups_change(attr, old, new):
            names = self._de_groups_checkbox.labels
            self._de_selected = {names[i] for i in new}

        def _on_clusters_change(attr, old, new):
            names = self._de_clusters_checkbox.labels
            self._de_selected_clusters = {names[i] for i in new}

        self._de_groups_checkbox.on_change("active", _on_groups_change)
        self._de_clusters_checkbox.on_change("active", _on_clusters_change)

        de_mode = RadioButtonGroup(labels=["vs Rest", "Between selected"], active=0)
        run_de_btn = Button(label="▶ Run DE", button_type="primary", width=120)
        level_select = self._de_level_select()
        self._de_logger = None
        logger = self._get_de_logger(doc)

        def on_run_de():
            if logger.running:
                return
            chosen = [
                (source, name)
                for source in self.DE_SOURCES
                for name in self._de_options(source)
                if name in self._de_chosen(source)
            ]
            if not chosen:
                logger.error("Select at least one group or cluster to compare.")
                return
            if de_mode.active == 1 and len(chosen) < 2:
                logger.error('"Between selected" needs at least 2 selected groups or clusters.')
                return

            names, masks, colors = [], {}, {}
            for source in self.DE_SOURCES:
                picked = [name for s, name in chosen if s == source]
                if picked:
                    names += picked
                    masks.update(self._de_masks(source, picked))
                    colors.update(self._de_colors(source))
            sources = {source for source, _ in chosen}
            self._run_de_thread(
                names,
                de_mode.active,
                doc,
                masks=masks,
                colors=colors,
                source=sources.pop() if len(sources) == 1 else "mixed",
                factor=self._de_factor_for(level_select),
            )

        run_de_btn.on_click(on_run_de)

        title = lambda text: Div(text=f"<span style='color:#888;font-size:11px'>{text}</span>")
        return [
            Div(text="<b style='color:#888888'>Compare</b>"),
            title("Groups"),
            self._de_groups_hint,
            self._de_groups_checkbox,
            title("Clusters"),
            self._de_clusters_hint,
            self._de_clusters_checkbox,
            Div(text="<b style='color:#888888'>Against</b>"),
            de_mode,
            *([level_select] if level_select is not None else []),
            run_de_btn,
            logger.div,
        ]

    def _clear_de_if_cluster_deleted(self):
        """Drop DE results when a referenced group/cluster is removed."""
        de_clusters = getattr(self, "_de_clusters", None)
        if not de_clusters or getattr(self, "_de_results_source", "groups") not in ("groups", "mixed"):
            return
        clusters = set(self._de_options("clusters"))
        if all(c in self._de_options("groups") or c in clusters for c in de_clusters):
            return
        self._clear_de_results()


__all__ = ["DEAnalysis", "DEManager"]
