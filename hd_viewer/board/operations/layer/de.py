
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from .... import settings
from ....startup.utils.progress import report

@dataclass
class DEResult:
    """Statistical results and the aligned expression/membership data used in testing.

    Expression is log1p-transformed after filtering genes. Its matrix is shared
    with the analysis cache and should be treated as read-only by consumers.
    """

    results: dict[str, pd.DataFrame]
    expression: ad.AnnData
    masks: dict[str, np.ndarray]
    names: tuple[str, ...]
    mode: int
    source: str
    factor: int
    n_overlap: int


class DEAnalysis:
    """Differential expression over explicit named observation masks.

    Instantiate with a sample exposing ``adata``; aggregated analysis additionally
    uses the sample's bin aggregation API. ``run`` is synchronous and can be used
    without GroupsManager, a component, a Bokeh document, or output files.
    """

    def __init__(self, sample):
        self.sample = sample
        self.prepared_data = {}

    def run(self, masks, mode=0, *, source="groups", factor=1):
        """Return a DEResult for each named set versus rest (0) or selected sets (1).

        Masks align with the sample's native observations. With factor > 1, counts
        are summed into bin aggregates and membership requires > 50% of their bins.
        Inputs are copied so preprocessing never modifies the sample's counts.
        """
        if mode not in (0, 1):
            raise ValueError("DE mode must be 0 (vs rest) or 1 (between selected).")
        if factor < 1 or int(factor) != factor:
            raise ValueError("Aggregation factor must be a positive integer.")
        if not masks:
            raise ValueError("Select at least one set to compare.")
        masks = {name: np.asarray(mask, dtype=bool).copy() for name, mask in masks.items()}
        if any(mask.shape != (self.sample.adata.n_obs,) for mask in masks.values()):
            raise ValueError("Each DE mask must have one value per sample observation.")
        active = list(masks)
        factor = int(factor)
        unit, units = {"groups": ("group", "groups"), "clusters": ("cluster", "clusters")}.get(
            source, ("set", "groups/clusters"))
        obs_word = "observations" if factor <= 1 else f"{factor} × {factor} bins"
        if factor > 1:
            level = self.sample._get_aggregated_sample(factor).adata
            level_masks = self.masks_at_level(masks, factor, level.obs_names)
        else:
            level, level_masks = self.sample.adata, masks

        cache = self.prepared_data
        if factor not in cache:
            report(f"Preparing expression data ({level.n_obs} {obs_word})")
            _a = level.copy()
            sc.pp.filter_genes(_a, min_cells=settings.MIN_GENE_CELLS)
            sc.pp.log1p(_a)
            cache[factor] = _a

        # filter_genes/log1p keep the observations: masks still line up
        adata_copy = cache[factor]
        results = {}
        tested_data = adata_copy
        n_overlap = 0

        if mode == 0:
            for cname in active:
                if cname not in level_masks:
                    continue
                in_cluster = pd.Series(level_masks[cname], index=adata_copy.obs.index)
                if (~in_cluster).sum() == 0 or in_cluster.sum() < settings.DE_MIN_OBSERVATIONS:
                    report(f"Skipping {cname}: fewer than {settings.DE_MIN_OBSERVATIONS} {obs_word}")
                    continue
                report(f"Testing {cname} vs rest")
                adata_copy.obs['_de_group'] = in_cluster.map(
                    {True: cname, False: 'Rest'}
                ).astype(str)
                try:
                    sc.tl.rank_genes_groups(
                        adata_copy, groupby='_de_group',
                        groups=[cname], reference='Rest',
                        method=settings.DE_METHOD,
                    )
                    results[cname] = sc.get.rank_genes_groups_df(adata_copy, group=cname)
                except Exception as exc:
                    report(f"{cname} failed: {exc}")
        else:
            # Label each observation with the selected set it belongs to.
            # Sets can overlap (groups with each other, a group with a cluster):
            # observations in several of them are left out, so each set is
            # compared on what only it contains. Independent of what the map shows.
            labels = pd.Series(None, index=adata_copy.obs.index, dtype=object)
            n_member = np.zeros(len(labels), dtype=int)
            for cname in active:
                if cname in level_masks:
                    n_member += level_masks[cname].astype(bool)
                    labels[level_masks[cname].astype(bool)] = cname
            overlap = n_member > 1
            labels[overlap] = None
            n_overlap = int(overlap.sum())
            if n_overlap:
                report(f"Left out {n_overlap} {obs_word} in more than one selected {unit}")
            in_sel = labels.notna()
            sub = adata_copy[in_sel.values].copy()
            sub.obs['_de_group'] = labels[in_sel].astype(str).values
            tested_data = sub

            counts = sub.obs['_de_group'].value_counts()
            usable = [c for c in active if counts.get(c, 0) >= settings.DE_MIN_OBSERVATIONS]
            if len(usable) < 2:
                raise ValueError(f"\"Between selected\" needs at least 2 {units} "
                                 f"with ≥{settings.DE_MIN_OBSERVATIONS} {obs_word} each"
                                 + (" outside their overlap." if n_overlap else "."))

            report(f"Testing {len(usable)} {units} against each other")
            sc.tl.rank_genes_groups(
                sub, groupby='_de_group',
                groups=usable, reference='rest',
                method=settings.DE_METHOD,
            )
            for cname in usable:
                results[cname] = sc.get.rank_genes_groups_df(sub, group=cname)

        if not results:
            raise ValueError(f"No DE results for the selected {units} "
                             f"(each needs ≥{settings.DE_MIN_OBSERVATIONS} {obs_word}).")

        comparison_masks = {
            name: (np.asarray(tested_data.obs['_de_group'].astype(str) == name, dtype=bool)
                   if mode == 1 else np.asarray(level_masks[name], dtype=bool).copy())
            for name in results
        }
        # Snapshot observation metadata without duplicating the large expression matrix.
        expression = ad.AnnData(X=tested_data.X,
                                obs=pd.DataFrame(index=tested_data.obs_names.copy()),
                                var=tested_data.var.copy())
        return DEResult(results=results, expression=expression, masks=comparison_masks,
                        names=tuple(active), mode=mode, source=source, factor=factor,
                        n_overlap=n_overlap)

    def masks_at_level(self, masks, factor, level_obs_names):
        """Native-bin memberships -> memberships of the aggregated cells of `factor`.

        A cell belongs to a set when most of its native bins do (> 50%). Clusters
        computed at the same level map exactly: all bins of a cell share its label.
        """
        cell_names = self.sample._aggregated_cell_names[factor]   # native obs -> cell
        valid = cell_names.notna().to_numpy()
        codes = pd.Categorical(cell_names[valid], categories=level_obs_names).codes
        n_bins = np.bincount(codes, minlength=len(level_obs_names))
        level_masks = {}
        for name, mask in masks.items():
            in_set = np.bincount(codes, weights=np.asarray(mask, dtype=bool)[valid].astype(float),
                                 minlength=len(level_obs_names))
            level_masks[name] = in_set > 0.5 * np.maximum(n_bins, 1)
        report(f"Mapped {len(masks)} sets to {len(level_obs_names)} aggregated cells "
               f"(member when > 50% of their bins are)")
        return level_masks


"""Native Bokeh views of the plot data prepared by DEPresenter."""

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Div, FixedTicker, HoverTool, Legend, LegendItem
from bokeh.plotting import figure


def build_de_plots(plot_data):
    """Build fresh models on the document thread; no image files or HTML embeds."""
    plots = []
    for data in plot_data:
        p = figure(title=f"{data.name} vs {data.rest_label} · top genes", height=340,
                   sizing_mode="stretch_width", min_width=320,
                   tools="pan,wheel_zoom,box_zoom,box_select,tap,reset,save",
                   active_scroll="wheel_zoom", background_fill_color="#14251f",
                   border_fill_color="#14251f", outline_line_color="#344d45")
        p.xaxis.ticker = FixedTicker(ticks=list(range(len(data.genes))))
        p.xaxis.major_label_overrides = dict(enumerate(data.genes))
        p.xaxis.major_label_orientation = .4
        p.yaxis.axis_label = "Log2 fold change" if data.kind == "dot" else "Expression · log1p(counts)"
        p.title.text_color = "#dddddd"
        p.axis.major_label_text_color = "#cccccc"
        p.axis.axis_label_text_color = "#cccccc"
        p.axis.axis_line_color = "#555555"
        p.grid.grid_line_color = "#333333"
        p.grid.grid_line_alpha = .4
        p.toolbar.logo = None
        hover_renderers, legend_items = [], []
        for group in dict.fromkeys(row["group"] for row in data.rows):
            rows = [row for row in data.rows if row["group"] == group]
            source = ColumnDataSource({key: [row[key] for row in rows] for key in rows[0]})
            if data.kind == "dot":
                dots = p.scatter("x", "lfc", source=source, size="size", color="color",
                                 alpha=.85, line_color="#eeeeee", selection_line_color="#ffffff",
                                 nonselection_alpha=.2, name="de_dots")
                hover_renderers.append(dots)
                legend_items.append(LegendItem(label=group, renderers=[dots]))
            else:
                violin = p.patches("xs", "ys", source=source, color="color", fill_alpha=.65,
                                   line_alpha=.9, nonselection_fill_alpha=.15, name="de_violins")
                quartiles = p.segment("x", "q1", "x", "q3", source=source,
                                      line_width=5, color="#dddddd")
                medians = p.scatter("x", "median", source=source, size=6, color="#ffffff")
                hover_renderers.extend([violin, medians])
                legend_items.append(LegendItem(label=group, renderers=[violin, quartiles, medians]))
        tooltips = [("Gene", "@gene"), ("Group", "@group"), ("Comparison", "@comparison"),
                    ("Log2 fold change", "@lfc{0.000}"), ("Adjusted p-value", "@padj"),
                    ("Expressing", "@pct{0.0}%"), ("Observations", "@n{0,0}"),
                    ("Mean log1p expression", "@mean{0.000}")]
        if data.kind == "violin":
            tooltips += [("Values in violin", "@n_shown{0,0}"), ("Q1", "@q1{0.000}"),
                         ("Median", "@median{0.000}"), ("Q3", "@q3{0.000}")]
        p.add_tools(HoverTool(renderers=hover_renderers, tooltips=tooltips))
        p.add_layout(Legend(items=legend_items, click_policy="hide", label_text_color="#cccccc",
                            background_fill_color="#192c28", border_line_color="#344d45"), "below")
        plots.append(p)
    return column(*plots, sizing_mode="stretch_width", spacing=16) if plots else column(
        Div(text="No positive marker genes to plot."), sizing_mode="stretch_width")


"""Presentation of DE results: table HTML and data for native Bokeh figures."""

from scipy.stats import gaussian_kde

@dataclass
class DEPlotData:
    """Plot-ready statistics and density coordinates; contains no UI models."""

    kind: str
    name: str
    rest_label: str
    genes: list[str]
    rows: list[dict]


@dataclass
class DEView:
    """Prepared presentation; construct Bokeh models on the document thread."""

    header: str
    merged: pd.DataFrame
    table_html: str
    plot_data: list[DEPlotData]


class DEPresenter:
    """Choose display labels, colours, marker genes and distribution geometry."""

    def __init__(self, alias="Spot", sample_type=None, colors=None):
        self.alias = alias
        self.sample_type = sample_type or f"{alias.lower()}s"
        self.colors = dict(colors or {})

    def _set_color(self, name):
        return self.colors.get(name) or '#4ef1a1'

    def prepare(self, result: DEResult) -> DEView:
        """Prepare view data off the UI thread without changing the analysis result."""
        self.result = result
        units = {"groups": "groups", "clusters": "clusters"}.get(result.source, "groups/clusters")
        keys = list(result.results)
        self._de_header = "DE vs Rest" if result.mode == 0 else (
            f"DE: {keys[0]} vs {keys[1]}" if len(keys) == 2 else f"DE between the selected {units}")
        if result.factor > 1:
            self._de_header += f" · {result.factor} × {result.factor} bins"
        report("Building table and plots")
        table_html = self._build_multi_de_table(result.results)
        plot_data = self._build_plot_data(result.names, result.expression, result.results, result.mode)
        return DEView(self._de_header, self._de_merged, table_html, plot_data)

    # Rows drawn at a time; the search box and sorting look through every gene
    DE_TABLE_ROWS = settings.DE_TABLE_ROWS

    # Inline handlers: scripts inserted in a Div do not run, event attributes do. Both
    # call the table's shared renderer (script.de-render), compiled once per table.
    _DE_RUN = ('(b._r||(b._r=new Function("b",'
               'b.querySelector("script.de-render").textContent)))(b);')
    _DE_SEARCH_JS = 'var b=this.closest("div.de-table");b._q=this.value;' + _DE_RUN
    _DE_SORT_JS = ('var b=this.closest("div.de-table");var c=+this.getAttribute("data-col");'
                   'if(b._col===c){b._dir=-b._dir;}else{b._col=c;'
                   'b._dir=this.getAttribute("data-first")==="desc"?-1:1;}' + _DE_RUN)

    def _build_multi_de_table(self, results: dict) -> str:
        all_genes: set = set()
        for df in results.values():
            all_genes.update(df['names'].tolist())

        merged = pd.DataFrame({'Gene': sorted(all_genes)})
        for cname, df in results.items():
            df = df.set_index('names')
            merged[f'{cname}_lfc']  = merged['Gene'].map(df['logfoldchanges']).round(2)
            merged[f'{cname}_pval'] = merged['Gene'].map(df['pvals_adj'])

        first_df = results[list(results.keys())[0]].set_index('names')
        merged['_sort'] = merged['Gene'].map(first_df['scores']).fillna(0)
        merged = merged.sort_values('_sort', ascending=False).drop(columns='_sort')

        self._de_merged = merged
        self._de_results = results
        return self._render_de_table_html(merged, results)

    @staticmethod
    def _de_cell_values(lfc, pval):
        """(text, colour) of the L2FC and p-adj cells; mirrored by the table's JS."""
        try:
            lfc_cell = (f"{float(lfc):.2f}", '#4ef1a1' if float(lfc) >= 0 else '#f14e4e')
            if np.isnan(float(lfc)):
                raise ValueError
        except (TypeError, ValueError):
            lfc_cell = ('—', '#888')
        try:
            if np.isnan(float(pval)):
                raise ValueError
            pval_cell = (f"{float(pval):.2e}", '#fff' if float(pval) < settings.DE_PVALUE_CUTOFF else '#666')
        except (TypeError, ValueError):
            pval_cell = ('—', '#666')
        return lfc_cell, pval_cell

    def _render_de_table_html(self, merged, results) -> str:
        import html as html_lib
        import json

        n_total = len(merged)
        names = list(results.keys())
        cluster_headers = ""
        for cname in names:
            color = self._set_color(cname)
            cluster_headers += f"""
                <th colspan='2' style='padding:5px 10px;color:{color};
                    font-weight:normal;text-align:center;
                    border-bottom:2px solid {color}'>{html_lib.escape(cname)}</th>"""

        # Click to sort: data-col indexes the JSON rows, data-first the first direction
        def sortable(label, col, first, align):
            return (f"<th data-col='{col}' data-first='{first}' onclick='{sort_attr}' "
                    f"title='Sort by {label}' style='padding:3px 8px;color:#888;text-align:{align};"
                    f"cursor:pointer;user-select:none'>{label}<span class='de-arrow'></span></th>")

        sort_attr = html_lib.escape(self._DE_SORT_JS, quote=True)
        gene_header = sortable("Gene", 0, "asc", "left")
        sub_headers = "".join(
            sortable("L2FC", 1 + 2 * i, "desc", "right") + sortable("p-adj", 2 + 2 * i, "asc", "right")
            for i in range(len(names))
        )

        def cell_html(text, color):
            return f"<td style='padding:4px 8px;text-align:right;color:{color}'>{text}</td>"

        rows_html = ""
        for _, row in merged.head(self.DE_TABLE_ROWS).iterrows():
            cells = f"<td style='padding:4px 10px;color:#fff'>{html_lib.escape(row['Gene'])}</td>"
            for cname in names:
                for text, color in self._de_cell_values(row.get(f'{cname}_lfc'), row.get(f'{cname}_pval')):
                    cells += cell_html(text, color)
            rows_html += f"<tr style='border-bottom:1px solid #222'>{cells}</tr>"

        # Every gene, for the search box: [gene, lfc, padj, lfc, padj, ...]; NaN -> null.
        # Only the displayed precision is sent (L2FC is already rounded to 2 decimals).
        def num(v, kind):
            if v is None or pd.isna(v):
                return None
            return float(v) if kind == 'lfc' else float(f"{float(v):.3g}")
        data = [[row['Gene']] + [num(row.get(f'{c}_{k}'), k) for c in names for k in ('lfc', 'pval')]
                for _, row in merged.iterrows()]
        data_json = json.dumps({"rows": data, "limit": self.DE_TABLE_ROWS, "pvalue_cutoff": settings.DE_PVALUE_CUTOFF},
                               separators=(',', ':'), allow_nan=False).replace("</", "<\\/")

        note = (f"Showing {min(n_total, self.DE_TABLE_ROWS)} of {n_total} genes"
                + (" · search to find any gene" if n_total > self.DE_TABLE_ROWS else ""))

        search_js = self._DE_SEARCH_JS
        render_js = (
            "var box=b;"
            "var d=box._de||(box._de=JSON.parse(box.querySelector(\"script.de-data\").textContent));"
            "var q=(box._q||\"\").trim().toLowerCase();"
            "var hits=q?d.rows.filter(function(r){return r[0].toLowerCase().indexOf(q)>=0;}):d.rows.slice();"
            "var c=box._col;"
            "if(c!=null){var s=box._dir;hits.sort(function(x,y){var a=x[c],z=y[c];"
            "if(a===null)return z===null?0:1;if(z===null)return -1;"
            "return c===0?s*a.localeCompare(z):s*(a-z);});}"
            "box.querySelectorAll(\"th[data-col]\").forEach(function(th){"
            "th.querySelector(\".de-arrow\").textContent=+th.getAttribute(\"data-col\")===c?(box._dir>0?\" ▲\":\" ▼\"):\"\";});"
            "function esc(t){var e=document.createElement(\"span\");e.textContent=t;return e.innerHTML;}"
            "function td(t,c){return \"<td style=\\\"padding:4px 8px;text-align:right;color:\"+c+\"\\\">\"+t+\"</td>\";}"
            "function exp(v){return v.toExponential(2).replace(/e([+-])(\\d)$/,\"e$10$2\");}"
            "var html=hits.slice(0,d.limit).map(function(r){"
            "var h=\"<td style=\\\"padding:4px 10px;color:#fff\\\">\"+esc(r[0])+\"</td>\";"
            "for(var i=1;i<r.length;i+=2){var l=r[i],p=r[i+1];"
            "h+=l===null?td(\"—\",\"#888\"):td(l.toFixed(2),l>=0?\"#4ef1a1\":\"#f14e4e\");"
            "h+=p===null?td(\"—\",\"#666\"):td(exp(p),p<d.pvalue_cutoff?\"#fff\":\"#666\");}"
            "return \"<tr style=\\\"border-bottom:1px solid #222\\\">\"+h+\"</tr>\";}).join(\"\");"
            "box.querySelector(\"tbody\").innerHTML=html;"
            "var n=Math.min(hits.length,d.limit);"
            "box.querySelector(\".de-note\").textContent=q"
            "?(hits.length+\" of \"+d.rows.length+\" genes match\"+(hits.length>n?\" (first \"+n+\" shown)\":\"\"))"
            ":(\"Showing \"+n+\" of \"+d.rows.length+\" genes\"+(d.rows.length>n?\" · search to find any gene\":\"\"));"
        )

        render_script = render_js.replace("</", "<\\/")   # keep the script block closed

        return f"""
        <div class='de-table' style='font-size:12px'>
            <script type='application/json' class='de-data'>{data_json}</script>
            <script type='text/plain' class='de-render'>{render_script}</script>
            <h4 style='color:#fff;margin:0 0 8px 0'>{html_lib.escape(getattr(self, '_de_header', 'DE vs Rest'))}</h4>
            <input
                type='text'
                placeholder='Search gene...'
                oninput='{html_lib.escape(search_js, quote=True)}'
                style='width:100%;box-sizing:border-box;margin-bottom:8px;
                    padding:5px 8px;background:#1a1a2e;color:#fff;
                    border:1px solid #444;border-radius:4px;font-size:12px'
            />
            <div style='overflow-y:auto;max-height:500px'>
                <table style='border-collapse:collapse;width:100%'>
                    <thead style='position:sticky;top:0;background:#000'>
                        <tr>
                            <th></th>
                            {cluster_headers}
                        </tr>
                        <tr>{gene_header}{sub_headers}</tr>
                    </thead>
                    <tbody>{rows_html}</tbody>
                </table>
            </div>
            <p class='de-note' style='color:#666;font-size:11px;margin:4px 0 0 0'>{note}</p>
        </div>"""

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------

    def _build_plot_data(self, active, adata, results, mode):
        plots = []
        for name in active:
            if name not in results:
                continue
            ranked = results[name]
            positive = ranked[ranked.logfoldchanges > 0].sort_values("scores", ascending=False)
            significant = positive[positive.pvals_adj < settings.DE_PVALUE_CUTOFF]
            top = pd.concat([significant, positive[~positive.names.isin(significant.names)]]).head(settings.DE_PLOT_TOP_GENES)
            if top.empty:
                continue
            if mode == 1:
                membership = np.asarray(self.result.masks[name], dtype=bool)
                rest_label = " + ".join(n for n in active if n != name) or "Others"
            else:
                membership = np.asarray(self.result.masks[name], dtype=bool)
                rest_label = "Rest"
            genes = top.names.tolist()
            rows = []
            for index, gene in enumerate(genes):
                values = adata[:, [gene]].X
                values = values.toarray().ravel() if hasattr(values, 'toarray') else np.asarray(values).ravel()
                stats = top.iloc[index]
                for offset, (group, mask, color) in enumerate((
                    (name, membership, self._set_color(name)),
                    (rest_label, ~membership, '#888888'),
                )):
                    expression = values[mask]
                    expression = expression[np.isfinite(expression)]
                    n = len(expression)
                    pct = 100 * float(np.count_nonzero(expression > 0)) / n if n else 0.
                    row = dict(gene=gene, group=group, comparison=f"{name} vs {rest_label}",
                               lfc=float(stats.logfoldchanges), padj=float(stats.pvals_adj),
                               pct=pct, n=n, mean=float(expression.mean()) if n else 0.,
                               color=color, size=4 + 20 * np.sqrt(pct / 100),
                               x=index + (-.18 if offset == 0 else .18))
                    if self.sample_type != 'spots':
                        # Preserve the existing positive-expression violin convention.
                        positive_values = expression[expression > 0]
                        shown = positive_values if len(positive_values) >= 3 else expression
                        row.update(n_shown=len(shown), q1=0., median=0., q3=0., xs=[], ys=[])
                        if len(shown):
                            q1, median, q3 = np.quantile(shown, [.25, .5, .75])
                            row.update(q1=float(q1), median=float(median), q3=float(q3))
                            row['xs'], row['ys'] = self._violin_shape(shown, row['x'])
                    rows.append(row)
            plots.append(DEPlotData('dot' if self.sample_type == 'spots' else 'violin',
                                    name, rest_label, genes, rows))
        return plots

    @staticmethod
    def _violin_shape(values, center):
        low, high = float(values.min()), float(values.max())
        if high <= low or len(values) < 2:
            delta = max(.02, abs(low) * .02)
            return ([center - .14, center + .14, center + .14, center - .14],
                    [max(0., low - delta), max(0., low - delta), high + delta, high + delta])
        grid = np.linspace(low, high, settings.DE_VIOLIN_GRID_POINTS)
        # Only the density approximation is sampled. Hover statistics use all values.
        sample = values if len(values) <= settings.DE_KDE_MAX_SAMPLES else np.random.default_rng(settings.DE_RANDOM_SEED).choice(values, settings.DE_KDE_MAX_SAMPLES, replace=False)
        try:
            density = gaussian_kde(sample)(grid)
        except np.linalg.LinAlgError:
            density = np.ones_like(grid)
        width = settings.DE_VIOLIN_HALF_WIDTH * density / max(float(density.max()), np.finfo(float).eps)
        return (np.concatenate([center - width, (center + width)[::-1]]).tolist(),
                np.concatenate([grid, grid[::-1]]).tolist())
