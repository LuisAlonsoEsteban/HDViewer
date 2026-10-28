"""Presentation of DE results: table HTML and data for native Bokeh figures."""

from .... import settings

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from .analysis import DEResult
from src.utils.progress import report


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
