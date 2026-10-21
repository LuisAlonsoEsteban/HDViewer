import threading
import pathlib
import uuid

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .loading import LoadingLogger
from src.utils.progress import report

GROUPS_COL = "custom_clusters"

# Filesystem path to bokeh_server/static/data — resolved relative to this file
# so it works on any machine without hardcoding.
# groups_de.py lives at:  bokeh_server/components/utils/groups_de.py
#   .parent → utils/   .parent → components/   .parent → bokeh_server/
_STATIC_DATA_DIR = pathlib.Path(__file__).resolve().parent.parent.parent / "static" / "data"
_STATIC_DATA_DIR.mkdir(parents=True, exist_ok=True)  # ignorado en git: crear si falta


class DEMixin:
    """DE analysis, HTML rendering, and plot generation for GroupsManager.

    Works on named sets of observations (custom groups or clusters), each given
    as a boolean mask aligned with the sample's observations.
    """

    def _set_color(self, name):
        lookup = getattr(self, '_de_color_lookup', None) or {}
        return lookup.get(name) or self._color_map.get(name, '#4ef1a1')

    # ------------------------------------------------------------------
    # DE thread
    # ------------------------------------------------------------------

    def _run_de_thread(self, active: list, mode: int, doc, masks=None, colors=None,
                       source="groups", factor=1) -> None:
        """mode 0: each set vs the rest; mode 1: between the selected sets.

        factor > 1 (bins): test the aggregated cells of that level (summed counts)
        instead of the native bins; see `_masks_at_level` for set membership.
        """
        logger = self._get_de_logger(doc)
        if logger.running:
            return
        # Words for messages: one source, or a mix of groups and clusters
        unit, units = {"groups": ("group", "groups"), "clusters": ("cluster", "clusters")}.get(
            source, ("set", "groups/clusters"))
        if masks is None:
            masks = self._de_masks("groups", active)
        self._de_masks_used = masks
        self._de_color_lookup = colors or {}
        self._de_results_source = source

        alias = self.component.alias.lower()
        # How the tested observations are called in messages, e.g. "bins" or "4 × 4 bins"
        obs_word = f"{alias}s" if factor <= 1 else f"{factor} × {factor} {alias}s"
        self._de_factor = factor
        # Conjuntos de este DE: si se borra alguno, se ocultan los resultados.
        self._de_clusters = list(active)

        def work():
            if factor > 1:
                level = self.sample._get_aggregated_sample(factor).adata
                level_masks = self._masks_at_level(masks, factor, level.obs_names)
            else:
                level, level_masks = self.sample.adata, masks
            self._de_masks_used = level_masks

            cache = self.__dict__.setdefault('_de_adata_by_factor', {})
            if factor not in cache:
                report(f"Preparing expression data ({level.n_obs} {obs_word})")
                _a = level.copy()
                sc.pp.filter_genes(_a, min_cells=3)
                sc.pp.log1p(_a)
                cache[factor] = _a

            # filter_genes/log1p keep the observations: masks still line up
            adata_copy = cache[factor]
            results = {}
            plot_adata = adata_copy
            n_overlap = 0

            if mode == 0:
                for cname in active:
                    if cname not in level_masks:
                        continue
                    in_cluster = pd.Series(level_masks[cname], index=adata_copy.obs.index)
                    if (~in_cluster).sum() == 0 or in_cluster.sum() < 3:
                        report(f"Skipping {cname}: fewer than 3 {obs_word}")
                        continue
                    report(f"Testing {cname} vs rest")
                    adata_copy.obs['_de_group'] = in_cluster.map(
                        {True: cname, False: 'Rest'}
                    ).astype(str)
                    try:
                        sc.tl.rank_genes_groups(
                            adata_copy, groupby='_de_group',
                            groups=[cname], reference='Rest',
                            method='wilcoxon',
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
                plot_adata = sub

                counts = sub.obs['_de_group'].value_counts()
                usable = [c for c in active if counts.get(c, 0) >= 3]
                if len(usable) < 2:
                    raise ValueError(f"\"Between selected\" needs at least 2 {units} "
                                     f"with ≥3 {obs_word} each"
                                     + (" outside their overlap." if n_overlap else "."))

                report(f"Testing {len(usable)} {units} against each other")
                sc.tl.rank_genes_groups(
                    sub, groupby='_de_group',
                    groups=usable, reference='rest',
                    method='wilcoxon',
                )
                for cname in usable:
                    results[cname] = sc.get.rank_genes_groups_df(sub, group=cname)

            if not results:
                raise ValueError(f"No DE results for the selected {units} "
                                 f"(each needs ≥3 {obs_word}).")

            _keys = list(results.keys())
            if mode != 0:
                self._de_header = (
                    f"DE: {_keys[0]} vs {_keys[1]}" if len(_keys) == 2
                    else f"DE between the selected {units}"
                )
            else:
                self._de_header = "DE vs Rest"
            if factor > 1:
                self._de_header += f" · {factor} × {factor} bins"

            report("Building table and plots")
            table_html = self._build_multi_de_table(results)
            plots_html = self._build_multi_plots_html(active, plot_adata, results, mode)
            self._last_plots_html = plots_html
            return len(results), n_overlap, table_html, plots_html

        def on_done(output):
            n_sets, n_overlap, table_html, plots_html = output
            if hasattr(self.component, '_figure_de_div'):
                self.component._figure_de_div.text = f"""
                <div style='display:grid;grid-template-columns:380px 1fr;gap:15px;width:100%'>
                    <div style='overflow-y:auto;max-height:600px'>{table_html}</div>
                    <div style='overflow-y:auto;max-height:600px'>{plots_html}</div>
                </div>"""
            left_out = f", {n_overlap} overlapping {obs_word} left out" if n_overlap else ""
            logger.success(f"DE done: {n_sets} {unit if n_sets == 1 else units}{left_out}, "
                           f"results below the map")

        logger.run(f"DE ({self._de_header_hint(mode, units)})", work,
                   on_done=on_done, on_error=lambda exc: logger.error(exc))

    def _masks_at_level(self, masks, factor, level_obs_names):
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

    def _de_header_hint(self, mode, units):
        return f"{units} vs rest" if mode == 0 else f"between selected {units}"

    def _get_de_logger(self, doc):
        """The DE tab's status line; created on demand (e.g. without the tab in tests)."""
        if getattr(self, '_de_logger', None) is None:
            self._de_logger = LoadingLogger(doc, width=270, font_size=11)
            self._de_status_div = self._de_logger.div
        return self._de_logger

    # ------------------------------------------------------------------
    # DE table
    # ------------------------------------------------------------------

    # Rows drawn at a time; the search box and sorting look through every gene
    DE_TABLE_ROWS = 200

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
            pval_cell = (f"{float(pval):.2e}", '#fff' if float(pval) < 0.05 else '#666')
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
        data_json = json.dumps({"rows": data, "limit": self.DE_TABLE_ROWS},
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
            "h+=p===null?td(\"—\",\"#666\"):td(exp(p),p<0.05?\"#fff\":\"#666\");}"
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

    def _build_multi_plots_html(
        self, active_clusters: list, adata, results: dict, mode: int = 0
    ) -> str:
        sample_type = getattr(self.component, 'sample_type', 'spots')
        if sample_type is None:
            raise ValueError(f"Component has no sample_type: {self.component}")

        static_dir = _STATIC_DATA_DIR

        def fig_to_url(fig):
            fname = f"de_plot_{uuid.uuid4().hex}.png"
            fig.savefig(str(static_dir / fname), format='png', dpi=100,
                        bbox_inches='tight', facecolor='#000000')
            plt.close(fig)
            return f"/bokeh_server/static/data/{fname}?t={uuid.uuid4().hex[:6]}"

        obs     = adata.obs.copy()
        unit_plural = f"{self.component.alias.lower()}s"   # spots / bins / cells
        n_genes = 5
        var_idx = {g: i for i, g in enumerate(adata.var_names)}

        GENE_COLORS = [
            '#4ef1a1', '#f1a14e', '#4e9af1', '#f14e4e', '#c44ef1',
            '#f1e44e', '#4ef1e4', '#f14ebc', '#8af14e', '#f17a4e',
        ]

        # Top genes per cluster
        cluster_genes: dict = {}
        for cname in active_clusters:
            if cname not in results:
                continue
            result = results[cname]
            sig = result[
                (result['pvals_adj'] < 0.05) & (result['logfoldchanges'] > 0)
            ].sort_values('scores', ascending=False)
            top = sig.head(n_genes)['names'].tolist()
            if len(top) < n_genes:
                rest = result[
                    (~result['names'].isin(top)) & (result['logfoldchanges'] > 0)
                ].sort_values('scores', ascending=False)
                top += rest.head(n_genes - len(top))['names'].tolist()
            cluster_genes[cname] = top

        violins_html = ""

        for cname in active_clusters:
            if cname not in cluster_genes or not cluster_genes[cname]:
                continue

            if mode != 0 and '_de_group' in obs.columns:
                in_cluster = obs['_de_group'].astype(str) == cname
                _others    = [c for c in active_clusters if c != cname]
                rest_label = " + ".join(_others) if _others else "Others"
            else:
                masks = getattr(self, '_de_masks_used', {})
                if cname not in masks:
                    continue
                in_cluster = pd.Series(masks[cname], index=obs.index)
                rest_label = "Rest"

            cluster_color = self._set_color(cname)
            top_genes     = cluster_genes[cname]
            n             = len(top_genes)

            if sample_type == 'spots':
                dot_pct_in, dot_pct_out, valid_genes = [], [], []

                for gene in top_genes:
                    if gene not in var_idx:
                        continue
                    expr = adata.X[:, var_idx[gene]]
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    in_vals  = np.log1p(expr[in_cluster.values])
                    out_vals = np.log1p(expr[~in_cluster.values])
                    dot_pct_in.append(np.mean(in_vals > 0) * 100)
                    dot_pct_out.append(np.mean(out_vals > 0) * 100)
                    valid_genes.append(gene)

                if not valid_genes:
                    continue

                lfc_map    = dict(zip(results[cname]['names'], results[cname]['logfoldchanges']))
                lfc_values = [lfc_map.get(g, 0) for g in valid_genes]
                order      = np.argsort(lfc_values)[::-1]
                dot_pct_in  = [dot_pct_in[i]  for i in order]
                dot_pct_out = [dot_pct_out[i] for i in order]
                valid_genes = [valid_genes[i]  for i in order]
                lfc_values  = [lfc_values[i]   for i in order]
                n_valid     = len(valid_genes)
                gene_colors = [GENE_COLORS[i % len(GENE_COLORS)] for i in range(n_valid)]

                fig, ax = plt.subplots(figsize=(max(7, 0.9 * n_valid + 3), 5))
                fig.patch.set_facecolor('#000000')
                ax.set_facecolor('#0a0a0a')

                y_pos    = np.arange(n_valid)
                max_size = 400

                for i, x in enumerate(y_pos):
                    bg = '#141414' if i % 2 == 0 else '#0a0a0a'
                    ax.axvspan(x - 0.5, x + 0.5, color=bg, zorder=0)
                ax.axhline(0, color='#333333', linewidth=0.8, zorder=1)

                for i, (gene, gc) in enumerate(zip(valid_genes, gene_colors)):
                    x     = y_pos[i]
                    s_in  = max_size * (dot_pct_in[i] / 100)
                    s_out = max_size * (dot_pct_out[i] / 100)
                    if lfc_values[i] > 0 or dot_pct_out[i] > 0:
                        ax.scatter(x - 0.18, lfc_values[i], s=max(s_out, 15), c='#444444',
                                   alpha=0.7, linewidths=0.8, edgecolors='#888888', zorder=3)
                        ax.plot([x - 0.18, x + 0.18], [lfc_values[i], lfc_values[i]],
                                color='#333333', linewidth=0.6, zorder=2)
                    ax.scatter(x + 0.18, lfc_values[i], s=max(s_in, 15), c=gc,
                               alpha=0.9, linewidths=0.8, edgecolors='white', zorder=4)

                ax.set_xticks(y_pos)
                ax.set_xticklabels(valid_genes, fontsize=12, rotation=30, ha='right')
                for tick, gc in zip(ax.get_xticklabels(), gene_colors):
                    tick.set_color(gc)
                ax.set_ylabel('Log2 Fold Change', color='white', fontsize=11)
                ax.set_xlabel('')
                ax.tick_params(colors='white', labelsize=7)
                plt.setp(ax.get_xticklabels(), rotation=20, ha='right')
                ax.set_ylim(bottom=-0.05)
                ax.set_xlim(-0.6, n_valid - 0.4)
                for spine in ['bottom', 'left']:
                    ax.spines[spine].set_color('#333')
                for spine in ['top', 'right']:
                    ax.spines[spine].set_visible(False)
                for label in ax.get_yticklabels():
                    label.set_color('white')

                from matplotlib.lines import Line2D
                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#aaaaaa',
                           markersize=8, label=f'{cname} (colored)', linestyle='None'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#444444',
                           markersize=6, label=rest_label, linestyle='None',
                           markeredgecolor='#888888'),
                ]
                leg1 = ax.legend(
                    handles=legend_elements, loc='upper left',
                    bbox_to_anchor=(0.9, 0.99), fontsize=10, framealpha=0.3,
                    labelcolor='white', facecolor='#111111', edgecolor='#444444',
                    title='Group', title_fontsize=10, borderaxespad=0,
                )
                plt.setp(leg1.get_title(), color='#aaaaaa')
                ax.add_artist(leg1)

                pct_levels   = [25, 50, 100]
                size_handles = [
                    plt.scatter([], [], s=max_size * p / 100, c='#888888', alpha=0.7,
                                edgecolors='white', linewidths=0.5, label=f'{p}%')
                    for p in pct_levels
                ]
                leg2 = ax.legend(
                    handles=size_handles, title='% expressing', title_fontsize=10,
                    fontsize=10, loc='upper left', bbox_to_anchor=(1.02, 0.5),
                    framealpha=0.3, labelcolor='white', facecolor='#111111',
                    edgecolor='#444444', handletextpad=0.3, borderpad=0.6,
                    borderaxespad=0,
                )
                plt.setp(leg2.get_title(), color='#aaaaaa')
                plt.tight_layout(pad=1.2)

                img_url = fig_to_url(fig)
                violins_html += f"""
                <div style='margin-bottom:16px'>
                    <p style='color:#aaa;font-size:11px;margin:0 0 4px 0'>
                        <span style='color:{cluster_color}'>●</span>
                        <b style='color:{cluster_color}'>{cname}</b> — top {n_valid} genes
                    </p>
                    <div style='overflow-x:auto;width:100%'>
                        <img src='{img_url}' style='width:500px;max-width:none;display:block'/>
                    </div>
                    <p style='color:#888;font-size:10px;margin:4px 0 0 0'>
                        X axis: top genes of {cname}. Y axis: log2 fold change ({cname} vs {rest_label}).
                        Dot size: % of {unit_plural} expressing the gene (colour = {cname}, grey = {rest_label}).
                    </p>
                </div>"""

            else:
                import seaborn as sns

                fig_w = max(2.5 * n, 5)
                fig, axes = plt.subplots(1, n, figsize=(fig_w, 3.5))
                if n == 1:
                    axes = [axes]
                fig.patch.set_facecolor('#000000')

                for ax, gene in zip(axes, top_genes):
                    ax.set_facecolor('#111111')
                    if gene not in var_idx:
                        continue
                    expr = adata.X[:, var_idx[gene]]
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    # adata es el del DE (_de_adata), que ya tiene log1p aplicado:
                    # no volver a transformar.
                    in_vals  = expr[in_cluster.values]
                    out_vals = expr[~in_cluster.values]

                    in_nz  = in_vals[in_vals > 0]   if (in_vals > 0).sum()  >= 3 else in_vals
                    out_nz = out_vals[out_vals > 0]  if (out_vals > 0).sum() >= 3 else out_vals

                    df_v = pd.DataFrame({
                        'Expression': np.concatenate([in_nz, out_nz]),
                        'Group':      [cname] * len(in_nz) + [rest_label] * len(out_nz),
                    })
                    sns.violinplot(
                        data=df_v, x='Group', y='Expression', ax=ax,
                        cut=0, inner='quartile', hue='Group', legend=False,
                        palette={cname: cluster_color, rest_label: '#555555'},
                        order=[cname, rest_label],
                    )

                    all_vals = np.concatenate([in_nz, out_nz])
                    ymax = np.nanpercentile(all_vals, 95) if len(all_vals) > 0 else 1.0
                    if np.isnan(ymax) or ymax <= 0:
                        ymax = 1.0
                    ax.set_ylim(0, ymax * 1.4)
                    ax.set_title(gene, color='white', fontsize=9)
                    ax.set_xlabel('')
                    ax.set_ylabel(
                        'Expression, log1p(counts)' if ax is axes[0] else '',
                        color='#aaaaaa', fontsize=9,
                    )
                    ax.tick_params(colors='#aaaaaa', labelsize=10)
                    for spine in ['bottom', 'left']:
                        ax.spines[spine].set_color('#444')
                    for spine in ['top', 'right']:
                        ax.spines[spine].set_visible(False)
                    for label in ax.get_xticklabels() + ax.get_yticklabels():
                        label.set_color('#aaaaaa')

                plt.tight_layout()
                img_url = fig_to_url(fig)
                img_w   = max(250 * n, 500)
                violins_html += f"""
                <div style='margin-bottom:16px'>
                    <p style='color:#aaa;font-size:11px;margin:0 0 4px 0'>
                        <span style='color:{cluster_color}'>●</span>
                        <b style='color:{cluster_color}'>{cname}</b> — top {n} genes
                    </p>
                    <div style='overflow-x:auto;width:100%'>
                        <img src='{img_url}' style='width:{img_w}px;max-width:none;display:block'/>
                    </div>
                    <p style='color:#888;font-size:10px;margin:4px 0 0 0'>
                        X axis: {cname} vs {rest_label}. Y axis: gene expression,
                        log1p(counts). Only {unit_plural} with expression &gt; 0 are shown;
                        the inner lines are the median and quartiles.
                    </p>
                </div>"""

        if not violins_html:
            violins_html = "<p style='color:#888;font-size:12px'>No results to display.</p>"

        return violins_html
