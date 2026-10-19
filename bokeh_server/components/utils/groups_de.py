import threading
import pathlib
import uuid

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GROUPS_COL = "custom_clusters"

# Filesystem path to bokeh_server/static/data — resolved relative to this file
# so it works on any machine without hardcoding.
# groups_de.py lives at:  bokeh_server/components/utils/groups_de.py
#   .parent → utils/   .parent → components/   .parent → bokeh_server/
_STATIC_DATA_DIR = pathlib.Path(__file__).resolve().parent.parent.parent / "static" / "data"
_STATIC_DATA_DIR.mkdir(parents=True, exist_ok=True)  # ignorado en git: crear si falta


class DEMixin:
    """DE analysis, HTML rendering, and plot generation for GroupsManager."""

    # ------------------------------------------------------------------
    # DE thread
    # ------------------------------------------------------------------

    def _run_de_thread(self, active: list, mode: int, doc) -> None:
        def _set_de_message(msg: str, color: str = "#aaa") -> None:
            # Todos los avisos/estado del DE van al _de_status_div, debajo del botón
            # Run DE (mismo sitio que los mensajes de bloqueo de on_run_de). Se
            # actualiza vía next tick para que sea seguro también desde el hilo y no
            # dejar la UI colgada en "Running DE..." cuando no hay resultados.
            if hasattr(self, '_de_status_div'):
                def _u():
                    self._de_status_div.text = (
                        f"<p style='color:{color};font-size:11px'>{msg}</p>"
                    )
                doc.add_next_tick_callback(_u)

        _set_de_message("Running DE...")
        # Clusters de este DE: si se deselecciona alguno, se ocultan los resultados.
        self._de_clusters = list(active)

        def _run():
            if getattr(self, '_de_adata', None) is None:
                _a = self.sample.adata.copy()
                sc.pp.filter_genes(_a, min_cells=3)
                sc.pp.log1p(_a)
                self._de_adata = _a

            adata_copy = self._de_adata
            for _c in self.sample.adata.obs.columns:
                if _c.startswith('grp_') or _c == GROUPS_COL:
                    adata_copy.obs[_c] = self.sample.adata.obs[_c]

            results = {}
            plot_adata = adata_copy

            if mode == 0:
                for cname in active:
                    col_key = f"grp_{cname}"
                    if col_key not in adata_copy.obs.columns:
                        continue
                    in_cluster = adata_copy.obs[col_key].astype(bool)
                    if (~in_cluster).sum() == 0 or in_cluster.sum() < 3:
                        continue
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
                    except Exception as e:
                        pass
            else:
                # Label each spot with the first selected group it belongs to.
                # GROUPS_COL only holds the groups visible on the map, and DE
                # must not depend on what is shown.
                labels = pd.Series(None, index=adata_copy.obs.index, dtype=object)
                for cname in active:
                    col_key = f"grp_{cname}"
                    if col_key in adata_copy.obs.columns:
                        member = adata_copy.obs[col_key].astype(bool) & labels.isna()
                        labels[member] = cname
                in_sel = labels.notna()
                sub = adata_copy[in_sel.values].copy()
                sub.obs['_de_group'] = labels[in_sel].astype(str).values
                plot_adata = sub

                counts = sub.obs['_de_group'].value_counts()
                usable = [c for c in active if counts.get(c, 0) >= 3]
                if len(usable) < 2:
                    _set_de_message(
                        "DE entre clusters necesita al menos 2 clusters con "
                        "≥3 spots cada uno.",
                        color="red",
                    )
                    return

                sc.tl.rank_genes_groups(
                    sub, groupby='_de_group',
                    groups=usable, reference='rest',
                    method='wilcoxon',
                )
                for cname in usable:
                    results[cname] = sc.get.rank_genes_groups_df(sub, group=cname)

            if not results:
                _set_de_message(
                    "DE sin resultados para la selección actual.", color="red"
                )
                return

            _keys = list(results.keys())
            if mode != 0:
                self._de_header = (
                    f"DE: {_keys[0]} vs {_keys[1]}" if len(_keys) == 2
                    else "DE entre clusters seleccionados"
                )
            else:
                self._de_header = "DE vs Rest"

            table_html = self._build_multi_de_table(results)
            plots_html = self._build_multi_plots_html(active, plot_adata, results, mode)
            self._last_plots_html = plots_html

            def _update():
                if hasattr(self, '_de_status_div'):
                    self._de_status_div.text = ""
                if hasattr(self.component, '_figure_de_div'):
                    self.component._figure_de_div.text = f"""
                    <div style='display:grid;grid-template-columns:380px 1fr;gap:15px;width:100%'>
                        <div style='overflow-y:auto;max-height:600px'>{table_html}</div>
                        <div style='overflow-y:auto;max-height:600px'>{plots_html}</div>
                    </div>"""

            doc.add_next_tick_callback(_update)

        t = threading.Thread(target=_run, daemon=True)
        t.start()

    # ------------------------------------------------------------------
    # DE table
    # ------------------------------------------------------------------

    def _build_multi_de_table(self, results: dict, filter_text: str = "") -> str:
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

        if filter_text:
            merged = merged[merged['Gene'].str.lower().str.contains(filter_text.lower())]
        else:
            merged = merged.head(200)

        return self._render_de_table_html(merged, results)

    def _render_de_table_html(self, merged, results, n_total=None) -> str:
        n_shown = len(merged)
        if n_total is None:
            n_total = len(self._de_merged) if self._de_merged is not None else n_shown

        cluster_headers = ""
        for cname in results.keys():
            color = self._color_map.get(cname, '#4ef1a1')
            cluster_headers += f"""
                <th colspan='2' style='padding:5px 10px;color:{color};
                    font-weight:normal;text-align:center;
                    border-bottom:2px solid {color}'>{cname}</th>"""

        sub_headers = "".join([
            "<th style='padding:3px 8px;color:#888;text-align:right'>L2FC</th>"
            "<th style='padding:3px 8px;color:#888;text-align:right'>p-adj</th>"
            for _ in results
        ])

        rows_html = ""
        for _, row in merged.iterrows():
            cells = f"<td style='padding:4px 10px;color:#fff'>{row['Gene']}</td>"
            for cname in results.keys():
                lfc  = row.get(f'{cname}_lfc',  float('nan'))
                pval = row.get(f'{cname}_pval', float('nan'))
                try:
                    lfc_color = '#4ef1a1' if float(lfc) >= 0 else '#f14e4e'
                    lfc_str   = f"{float(lfc):.2f}"
                except Exception:
                    lfc_color = '#888'
                    lfc_str   = '—'
                try:
                    pval_str   = f"{float(pval):.2e}"
                    pval_color = '#fff' if float(pval) < 0.05 else '#666'
                except Exception:
                    pval_str   = '—'
                    pval_color = '#666'
                cells += (
                    f"<td style='padding:4px 8px;text-align:right;color:{lfc_color}'>{lfc_str}</td>"
                    f"<td style='padding:4px 8px;text-align:right;color:{pval_color}'>{pval_str}</td>"
                )
            rows_html += (
                f"<tr class='de-row' data-gene='{row['Gene'].lower()}'"
                f" style='border-bottom:1px solid #222'>{cells}</tr>"
            )

        return f"""
        <div style='font-size:12px'>
            <h4 style='color:#fff;margin:0 0 8px 0'>{getattr(self, '_de_header', 'DE vs Rest')}</h4>
            <input
                id='gene-search-input'
                type='text'
                placeholder='Search gene...'
                onkeyup='
                    var q = this.value.toLowerCase();
                    var rows = this.parentNode.querySelectorAll("tr.de-row");
                    for (var i=0; i<rows.length; i++) {{
                        var gene = rows[i].getAttribute("data-gene");
                        rows[i].style.display = gene.indexOf(q) >= 0 ? "" : "none";
                    }}
                '
                style='width:100%;box-sizing:border-box;margin-bottom:8px;
                    padding:5px 8px;background:#1a1a2e;color:#fff;
                    border:1px solid #444;border-radius:4px;font-size:12px'
            />
            <div style='overflow-y:auto;max-height:500px'>
                <table style='border-collapse:collapse;width:100%'>
                    <thead style='position:sticky;top:0;background:#000'>
                        <tr>
                            <th style='padding:5px 10px;color:#888;text-align:left'>Gene</th>
                            {cluster_headers}
                        </tr>
                        <tr><th></th>{sub_headers}</tr>
                    </thead>
                    <tbody>{rows_html}</tbody>
                </table>
            </div>
            <p style='color:#666;font-size:11px;margin:4px 0 0 0'>
                Showing {n_shown} of {n_total} genes
            </p>
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
                rest_label = " + ".join(_others) if _others else "Otros"
            else:
                col_key = f"grp_{cname}"
                if col_key not in obs.columns:
                    continue
                in_cluster = obs[col_key].astype(bool)
                rest_label = "Rest"

            cluster_color = self._color_map.get(cname, '#4ef1a1')
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
                        Eje X: top genes de {cname}. Eje Y: log2 fold change ({cname} vs {rest_label}).
                        Tamaño del punto: % de spots que expresan el gen (color = {cname}, gris = {rest_label}).
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
                        'Expresión log1p(conteos)' if ax is axes[0] else '',
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
                        Eje X: grupo ({cname} vs {rest_label}). Eje Y: expresión del gen,
                        log1p(conteos). Solo se muestran las observaciones con expresión &gt; 0;
                        las líneas internas son la mediana y los cuartiles.
                    </p>
                </div>"""

        if not violins_html:
            violins_html = "<p style='color:#888;font-size:12px'>No results to display.</p>"

        return violins_html
