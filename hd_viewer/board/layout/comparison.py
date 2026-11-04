
import html

import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    ColumnDataSource,
    DataTable,
    Div,
    HelpButton,
    HTMLTemplateFormatter,
    NumberFormatter,
    Select,
    Slider,
    TableColumn,
    TabPanel,
    Tabs,
    Tooltip,
)
from bokeh.models.dom import HTML
from bokeh.plotting import figure

from ... import settings
from ...settings.config import comparison_kwargs
from . import Panel
from ..utils.loading import LoadingLogger
from ...help.comparison import GENE_HELP, LOCAL_HELP, global_help


class ComparisonPanel(Panel):
    """Compare a low-resolution and a high-resolution sample.

    Per-geometry metrics are written to the low-res `adata.obs` as `cmp_<metric>` columns, so they
    can be shown on the map like any other metadata column.
    """

    OBS_PREFIX = settings.COMPARISON_OBS_PREFIX
    GLOBAL_HELP_HINT = "<i>&#9432; Click a metric to see what it captures.</i>"
    NONE = "(none)"

    GENE_COLUMNS = [
        ("spatial_pearson", "Spatial r"),
        ("capture_ratio", "Capture"),
        ("morans_i_sample1", "Moran low"),
        ("morans_i_sample2", "Moran high"),
        ("zero_fraction_sample1", "Zeros low"),
        ("zero_fraction_sample2", "Zeros high"),
    ]

    def __init__(self, doc):
        super().__init__(name="comparison_panel")
        self.doc = doc

        self.components = {}
        self.result = None
        # Whether the user's Sample 1 is the high-res one (result names are in the user's order,
        # help texts and the gene table use the canonical low-res = sample1 names)
        self.swapped = False
        self.sample1_component = None
        self.sample2_component = None

        self._build_controls()
        self._build_results()

        self.set_children([
            self.sample1_select,
            self.sample2_select,
            row(self.predicate_select, self.coverage_slider, sizing_mode="stretch_width"),
            row(self.sample1_cluster_select, self.sample2_cluster_select, sizing_mode="stretch_width"),
            row(self.sample1_celltype_select, self.sample2_celltype_select, sizing_mode="stretch_width"),
            self.run_button,
            self.status.div,
            self.result_tabs,
        ])

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _build_controls(self):
        self.sample1_select = Select(title="Sample 1", options=[])
        self.sample2_select = Select(title="Sample 2", options=[])
        self.sample1_select.on_change("value", lambda attr, old, new: self._refresh_key_options())
        self.sample2_select.on_change("value", lambda attr, old, new: self._refresh_key_options())

        self.predicate_select = Select(title="Match high-res geometries that are",
                                       options=[("contains", "inside"), ("intersects", "overlapping")],
                                       value=settings.COMPARISON_PREDICATE)
        self.coverage_slider = Slider(title="Min. coverage", start=0.0, end=1.0, step=settings.COMPARISON_COVERAGE_STEP, value=settings.COMPARISON_MIN_COVERAGE)

        self.sample1_cluster_select = Select(title="Sample 1 clusters", options=[self.NONE], value=self.NONE)
        self.sample2_cluster_select = Select(title="Sample 2 clusters", options=[self.NONE], value=self.NONE)
        self.sample1_celltype_select = Select(title="Sample 1 cell types", options=[self.NONE], value=self.NONE)
        self.sample2_celltype_select = Select(title="Sample 2 cell types", options=[self.NONE], value=self.NONE)

        self.run_button = Button(label="Compare samples", button_type="success", disabled=True)
        self.run_button.on_click(self._on_run)

        self.status = LoadingLogger(self.doc, text="Load at least two samples to compare them.")

    def _build_results(self):
        # Global metrics
        self.global_source = ColumnDataSource(data={"metric": [], "value": [], "help": []})
        # Clicking a row explains its metric below the table (the icon also has a hover tooltip)
        self.global_source.selected.on_change("indices", self._on_global_selection)
        self.global_help_div = Div(text=self.GLOBAL_HELP_HINT, sizing_mode="stretch_width",
                                   styles={"padding": "6px", "border-left": "3px solid #4e9af1"})
        metric_formatter = HTMLTemplateFormatter(
            template='<span title="<%= help %>" style="cursor:help;color:#4e9af1">&#9432;</span> <%= value %>')
        global_table = DataTable(
            source=self.global_source,
            columns=[TableColumn(field="metric", title="Metric", formatter=metric_formatter),
                     TableColumn(field="value", title="Value", formatter=NumberFormatter(format="0.000"))],
            sizing_mode="stretch_width", height=380, index_position=None,
        )

        # Local metrics: show one on the map, its histogram and a scatter of two metrics
        self.metric_select = Select(title="Local metric", options=[])
        self.metric_select.on_change("value", lambda attr, old, new: self._update_histogram())
        self.metric_help = Tooltip(content="Choose a metric.", position="right")
        self.metric_select.on_change("value", lambda attr, old, new: self._update_metric_help())
        metric_help_button = HelpButton(tooltip=self.metric_help, align="end")
        self.show_metric_button = Button(label="Show on map")
        self.show_metric_button.on_click(self._show_metric_on_map)

        self.hist_source = ColumnDataSource(data={"top": [], "left": [], "right": []})
        self.hist_plot = figure(height=200, sizing_mode="stretch_width", tools="pan,wheel_zoom,reset,save")
        self.hist_plot.quad(top="top", bottom=0, left="left", right="right", source=self.hist_source,
                            fill_color="#4e9af1", line_color="#1e1e2e", alpha=0.8)
        self.hist_plot.toolbar.logo = None

        self.scatter_x_select = Select(title="x", options=[])
        self.scatter_y_select = Select(title="y", options=[])
        self.scatter_x_select.on_change("value", lambda attr, old, new: self._update_scatter())
        self.scatter_y_select.on_change("value", lambda attr, old, new: self._update_scatter())

        self.scatter_source = ColumnDataSource(data={"x": [], "y": [], "obs_name": []})
        self.scatter_source.selected.on_change("indices", self._on_scatter_selection)
        self.scatter_plot = figure(height=250, sizing_mode="stretch_width",
                                   tools="pan,wheel_zoom,box_select,lasso_select,reset,save",
                                   tooltips=[("", "@obs_name"), ("x", "@x"), ("y", "@y")])
        self.scatter_plot.scatter(x="x", y="y", source=self.scatter_source, size=4, alpha=0.6, color="#f1a14e")
        self.scatter_plot.toolbar.logo = None

        local_tab = column(
            row(self.metric_select, metric_help_button, self.show_metric_button, sizing_mode="stretch_width"),
            self.hist_plot,
            row(self.scatter_x_select, self.scatter_y_select, sizing_mode="stretch_width"),
            Div(text="<i>Select points to highlight the geometries on the map.</i>"),
            self.scatter_plot,
            sizing_mode="stretch_width",
        )

        # Per-gene metrics: clicking a gene shows it on both samples
        self.gene_source = ColumnDataSource(data={"gene": [], **{col: [] for col, _ in self.GENE_COLUMNS}})
        self.gene_source.selected.on_change("indices", self._on_gene_selection)
        gene_table = DataTable(
            source=self.gene_source,
            columns=[TableColumn(field="gene", title="Gene")] + [
                TableColumn(field=col, title=title, formatter=NumberFormatter(format="0.000"))
                for col, title in self.GENE_COLUMNS
            ],
            sizing_mode="stretch_width", height=380, index_position=None, sortable=True,
        )
        gene_help = "<br>".join(f"<b>{title}</b>: {html.escape(GENE_HELP.get(col, ''))}"
                                for col, title in self.GENE_COLUMNS)
        gene_help_button = HelpButton(tooltip=Tooltip(content=HTML(gene_help), position="right"))
        genes_tab = column(row(Div(text="<i>Click a gene to show it on both samples.</i>"), gene_help_button),
                           gene_table, sizing_mode="stretch_width")

        self.result_tabs = Tabs(tabs=[
            TabPanel(child=column(global_table, self.global_help_div, sizing_mode="stretch_width"), title="Global"),
            TabPanel(child=local_tab, title="Local"),
            TabPanel(child=genes_tab, title="Genes"),
        ], visible=False, sizing_mode="stretch_width")

    # ------------------------------------------------------------------
    # Samples
    # ------------------------------------------------------------------

    def set_components(self, components):
        self.components = {c.name: c for c in components if c.sample.has_adata()}
        aliases = list(self.components.keys())

        # Default order: largest geometries are the low-res sample
        by_area = sorted(aliases, key=lambda a: -self.components[a].sample.gdf.geometry.area.mean())
        for select, default in ((self.sample1_select, by_area[:1]), (self.sample2_select, by_area[1:2])):
            select.options = aliases
            if select.value not in aliases and default:
                select.value = default[0]

        self.run_button.disabled = len(aliases) < 2
        if len(aliases) >= 2 and self.result is None:
            self.status.div.text = "Choose the samples and press <b>Compare samples</b>."
        self._refresh_key_options()

    def _refresh_key_options(self):
        # Options follow the user's order: compare_samples moves each key together with its sample
        sample1 = self.components.get(self.sample1_select.value)
        sample2 = self.components.get(self.sample2_select.value)
        self._set_options(self.sample1_cluster_select, self._label_columns(sample1))
        self._set_options(self.sample2_cluster_select, self._label_columns(sample2))
        self._set_options(self.sample1_celltype_select, self._celltype_options(sample1))
        self._set_options(self.sample2_celltype_select, self._celltype_options(sample2))

    def _set_options(self, select, options):
        select.options = [self.NONE] + options
        values = [o[0] if isinstance(o, tuple) else o for o in select.options]
        if select.value not in values:
            select.value = self.NONE

    def _celltype_options(self, component):
        """Deconvolution proportions (obsm DataFrames) first, then label columns of adata.obs."""
        proportions = self._proportion_keys(component)
        labels = [c for c in self._label_columns(component) if c not in proportions]
        return [(k, f"{k} (proportions)") for k in proportions] + [(c, c) for c in labels]

    def _label_columns(self, component):
        if component is None:
            return []
        obs = component.sample.adata.obs
        return [c for c in obs.columns
                if not c.startswith(self.OBS_PREFIX) and not pd.api.types.is_numeric_dtype(obs[c])]

    @staticmethod
    def _proportion_keys(component):
        if component is None:
            return []
        obsm = component.sample.adata.obsm
        return [k for k in obsm.keys() if isinstance(obsm[k], pd.DataFrame)]

    def _key(self, select):
        return None if select.value == self.NONE else select.value

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def _on_run(self):
        low = self.components.get(self.sample1_select.value)
        high = self.components.get(self.sample2_select.value)
        if low is None or high is None or low is high:
            self.status.set_message("Choose two different samples.", color="red")
            return

        kwargs = dict(
            **{k: v for k, v in comparison_kwargs().items() if k not in ("min_coverage", "predicate")},
            min_coverage=self.coverage_slider.value,
            predicate=self.predicate_select.value,
            sample1_cluster_key=self._key(self.sample1_cluster_select),
            sample2_cluster_key=self._key(self.sample2_cluster_select),
            sample1_celltype_key=self._key(self.sample1_celltype_select),
            sample2_celltype_key=self._key(self.sample2_celltype_select),
        )

        if self.status.running:
            return

        def _work():
            from ..operations.comparison import compare_samples
            return compare_samples(low.sample, high.sample, **kwargs)

        self.run_button.disabled = True
        self.status.run(
            f"Comparing {low.name} with {high.name}",
            _work,
            on_done=lambda result: self._on_result(low, high, result),
            on_error=self._on_error,
        )

    def _on_error(self, exc):
        self.status.error(exc)
        self.run_button.disabled = False

    def _on_result(self, sample1, sample2, result):
        self.run_button.disabled = False

        # Local metrics live on the low-res (reference) sample; names stay in the user's order
        self.swapped = result.get("reference_sample", 1) == 2
        low, high = (sample2, sample1) if self.swapped else (sample1, sample2)

        # Results of a previous run are overwritten (or dropped when nothing matches)
        shown = {}
        for component in {self.sample1_component, low} - {None}:
            shown[component] = self._clear_obs_columns(component)

        self.sample1_component, self.sample2_component = low, high
        self.result = result

        self.global_source.selected.indices = []
        self.global_source.data = self._flatten_global(result["global"])
        self.global_help_div.text = self.GLOBAL_HELP_HINT
        self.result_tabs.visible = True

        local = result["local"]
        if local is None:
            self.status.set_message(f"No {low.alias.lower()} passes the coverage filter: "
                                    "only tissue overlap metrics are available.", color="orange")
            for source in (self.hist_source, self.scatter_source, self.gene_source):
                source.data = {key: [] for key in source.data}

            self.metric_select.options = []
            low.refresh_completions()
            return

        metrics = self._write_obs_columns(low, local)
        numeric = [m for m in metrics if pd.api.types.is_numeric_dtype(local[m])]
        self.metric_select.options = numeric

        self.scatter_x_select.options = numeric
        self.scatter_y_select.options = numeric

        self.metric_select.value = "pearson_hvg" if "pearson_hvg" in numeric else numeric[0]

        self.scatter_x_select.value = "total_counts_sample1"
        self.scatter_y_select.value = "total_counts_sample2"

        self._update_histogram()
        self._update_scatter()

        genes = result["genes"].sort_values("spatial_pearson", ascending=False)
        self.gene_source.data = {
            "gene": genes.index.astype(str).tolist(),
            **{col: genes[self._name(col)].astype(float).tolist() for col, _ in self.GENE_COLUMNS},
        }

        low.refresh_completions()

        if shown.get(low):
            self._show_column(low, shown[low])
        self.status.success(
            f"Compared {len(local)} {low.alias.lower()}s with "
            f"{int(result['global'][self._name('n_sample2_matched')])} {high.alias.lower()}s. Local metrics are available "
            f"as <code>{self.OBS_PREFIX}*</code> columns of the {low.alias.lower()}s."
        )

    def _name(self, canonical):
        """Result name of a canonical (low-res = sample1) metric name, and vice versa."""
        from ..operations.comparison import oriented_name
        return oriented_name(canonical, self.swapped)

    def _flatten_global(self, global_metrics):
        metrics, values = [], []
        for key, value in global_metrics.items():
            items = value.items() if isinstance(value, dict) else [(None, value)]
            for sub, val in items:
                metrics.append(key if sub is None else f"{key} [{sub}]")
                values.append(float(val) if val is not None else np.nan)
        return {"metric": metrics, "value": values,
                "help": [html.escape(global_help(self._name(m)), quote=True) for m in metrics]}

    def _write_obs_columns(self, component, local):
        """Store the local metrics in the low-res adata.obs, NaN (not sequenced) where unmatched."""
        obs = component.sample.adata.obs
        obs_names = obs.index.astype(str)

        metrics = []
        for col in local.columns:
            # The low-res sample's own clusters are already one of its columns
            if col in ("obs_name", self._name("cluster_sample1")):
                continue
            values = pd.Series(local[col].to_numpy(), index=local["obs_name"].to_numpy())
            if not pd.api.types.is_numeric_dtype(values) or pd.api.types.is_bool_dtype(values):
                values = values.map(lambda v: str(v) if pd.notna(v) else np.nan)
            else:
                values = values.astype(float).replace([np.inf, -np.inf], np.nan)
            obs[self.OBS_PREFIX + col] = values.reindex(obs_names).to_numpy()
            metrics.append(col)
        return metrics

    def _clear_obs_columns(self, component):
        """Drop the metric columns of a previous run, returning the one on display (if any)."""
        shown = None
        autocomplete = component.autocomplete_input
        if autocomplete is not None and autocomplete.value.startswith(self.OBS_PREFIX):
            shown = autocomplete.value
            autocomplete.value = component.get_colnames()[0]

        cache = component.manager.cache
        cache.drop(columns=[c for c in cache.columns if str(c).startswith(self.OBS_PREFIX)], inplace=True)

        obs = component.sample.adata.obs
        obs.drop(columns=[c for c in obs.columns if c.startswith(self.OBS_PREFIX)], inplace=True)
        return shown


    # ------------------------------------------------------------------
    # Interaction
    # ------------------------------------------------------------------

    def _local_values(self, metric):
        local = self.result["local"] if self.result else None
        if local is None or metric not in local.columns:
            return None
        return pd.to_numeric(local[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)

    def _on_global_selection(self, attr, old, new):
        data = self.global_source.data
        if not new or new[0] >= len(data["metric"]):
            self.global_help_div.text = self.GLOBAL_HELP_HINT
            return
        i = new[0]
        self.global_help_div.text = f"<b>{data['metric'][i]}</b>: {data['help'][i] or 'No description.'}"

    def _update_metric_help(self):
        metric = self.metric_select.value
        self.metric_help.content = LOCAL_HELP.get(self._name(metric), "") if metric else "Choose a metric."

    def _update_histogram(self):
        values = self._local_values(self.metric_select.value)
        if values is None or values.notna().sum() == 0:
            self.hist_source.data = {"top": [], "left": [], "right": []}
            return
        hist, edges = np.histogram(values.dropna(), bins=30)
        self.hist_source.data = {"top": hist.tolist(), "left": edges[:-1].tolist(), "right": edges[1:].tolist()}
        self.hist_plot.xaxis.axis_label = self.metric_select.value
        self.hist_plot.yaxis.axis_label = f"N {self.sample1_component.alias.lower()}s"

    def _update_scatter(self):
        x = self._local_values(self.scatter_x_select.value)
        y = self._local_values(self.scatter_y_select.value)
        if x is None or y is None:
            self.scatter_source.data = {"x": [], "y": [], "obs_name": []}
            return
        ok = (x.notna() & y.notna()).to_numpy()
        self.scatter_source.selected.indices = []
        self.scatter_source.data = {
            "x": x[ok].tolist(),
            "y": y[ok].tolist(),
            "obs_name": self.result["local"]["obs_name"].to_numpy()[ok].tolist(),
        }
        self.scatter_plot.xaxis.axis_label = self.scatter_x_select.value
        self.scatter_plot.yaxis.axis_label = self.scatter_y_select.value

    def _show_column(self, component, colname):
        """Display `colname` on a sample through its Gene/Metadata input."""
        if component is None or component.autocomplete_input is None:
            return

        component.refresh_completions()

        if colname in component.get_colnames():
            component.autocomplete_input.value = colname

    def _show_metric_on_map(self):
        if self.metric_select.value:
            self._show_column(self.sample1_component, self.OBS_PREFIX + self.metric_select.value)

    def _on_scatter_selection(self, attr, old, new):
        if self.sample1_component is None:
            return
        obs_names = self.scatter_source.data["obs_name"]
        self.sample1_component.set_selected_obs_names([obs_names[i] for i in new if i < len(obs_names)])

    def _on_gene_selection(self, attr, old, new):
        if not new:
            return
        gene = self.gene_source.data["gene"][new[0]]
        self._show_column(self.sample1_component, gene)
        self._show_column(self.sample2_component, gene)
