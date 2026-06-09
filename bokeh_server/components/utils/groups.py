import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.models import (
    Button, Div, Select, TextInput,
    AutocompleteInput, CheckboxGroup,
)
from bokeh.layouts import column, row

from bokeh.models import Button, Div, Select

OPERATORS = [">", "<", ">=", "<=", "=="]

# Columna fija en adata.obs, igual que 'leiden'
GROUPS_COL = "custom_clusters"


class GroupsManager:
    """
    Manages creation of custom clusters based on gene/metadata conditions.
    Stores results in adata.obs[GROUPS_COL], same structure as leiden.

    Each spot gets the name of the first matching cluster, or None if unassigned.
    """

    def __init__(self, component):
        self.component = component
        self._clusters = {}       # cluster_name → [conditions]
        self._conditions = []     # conditions being built [(col, op, val, logic)]

        self._active_clusters = set()
        self._checkbox_group = None
        self._cluster_layout = None  # updated in-place like ClusteringManager
        self._color_map = {} 
        self._conditions_div = None
        self._status_div = None

    @property
    def sample(self):
        return self.component.manager.sample

    @property
    def is_clustered(self):
        return (
            GROUPS_COL in self.sample.adata.obs.columns
            and len(self._clusters) > 0
        )

    # ------------------------------------------------------------------
    # Public: get_widgets
    # ------------------------------------------------------------------

    def get_widgets(self, doc) -> list:

        # --- Condition builder ---
        colnames = self.sample.get_colnames()

        gene_input = AutocompleteInput(
            title="Gene / Metadata",
            value=colnames[0] if colnames else "",
            completions=colnames,
            min_characters=1,
            case_sensitive=False,
            restrict=False,
            width=160
        )

        operator_select = Select(
            title="Operator",
            options=OPERATORS,
            value=">",
            width=70
        )

        value_input = TextInput(
            title="Value",
            value="1",
            width=80
        )

        logic_select = Select(
            title="Logic",
            options=["AND", "OR"],
            value="AND",
            width=70
        )

        self._conditions_div = Div(
            text="<p style='color:#888;font-size:11px'>No conditions added yet.</p>",
            width=380
        )

        cluster_name_input = TextInput(
            value="Cluster A",
            title="Cluster name",
            width=200
        )

        add_condition_btn = Button(
            label="+ Add condition",
            button_type="success",
            width=150
        )
        clear_conditions_btn = Button(
            label="Clear",
            button_type="warning",
            width=80
        )
        create_cluster_btn = Button(
            label="✓ Create cluster",
            button_type="success",
            width=160
        )
        delete_cluster_btn = Button(
            label="🗑 Delete cluster",
            button_type="danger",
            width=160
        )

        self._status_div = Div(text="", width=380)

        # Cluster list with checkboxes — updated in-place like ClusteringManager
        if self._cluster_layout is None:
            self._cluster_layout = column(
                name="groups_cluster_controls",
                children=self._build_cluster_widgets(doc),
            )
        

        # ------------------------------------------------------------------
        # Render helpers
        # ------------------------------------------------------------------

        def _render_conditions():
            if not self._conditions:
                self._conditions_div.text = (
                    "<p style='color:#888;font-size:11px'>No conditions added yet.</p>"
                )
                return
            rows = ""
            for i, (col, op, val, logic) in enumerate(self._conditions):
                prefix = f"<b>{logic}</b> " if i > 0 else ""
                rows += (
                    f"<div style='margin:2px 0;color:#dddddd'>"
                    f"{prefix}<code>{col} {op} {val}</code></div>"
                )
            self._conditions_div.text = (
                f"<div style='padding:5px;background:#2a2a3e;border-radius:4px'>{rows}</div>"
            )

        # ------------------------------------------------------------------
        # Callbacks
        # ------------------------------------------------------------------

        def on_add_condition():
            col = gene_input.value.strip()
            op = operator_select.value
            val = value_input.value.strip()
            logic = logic_select.value if self._conditions else "AND"

            if not col or not val:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Gene and value are required.</p>"
                )
                return
            try:
                float(val)
            except ValueError:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Value must be numeric.</p>"
                )
                return

            # Show gene range info
            if col in self.sample.adata.var_names:
                col_idx = list(self.sample.adata.var_names).index(col)
                expr = self.sample.adata.X[:, col_idx]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                self._status_div.text = (
                    f"<p style='color:#4e9af1;font-size:11px'>"
                    f"{col}: max={float(expr.max()):.2f}, "
                    f"mean={float(expr.mean()):.3f}</p>"
                )
            else:
                self._status_div.text = ""

            self._conditions.append((col, op, val, logic))
            doc.add_next_tick_callback(_render_conditions)

        def on_clear_conditions():
            self._conditions.clear()
            doc.add_next_tick_callback(_render_conditions)

        def on_create_cluster():
            if not self._conditions:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Add at least one condition first.</p>"
                )
                return

            cname = cluster_name_input.value.strip() or "Cluster A"

            if cname in self._clusters:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>"
                    f"Cluster '{cname}' already exists.</p>"
                )
                return

            self._clusters[cname] = list(self._conditions)
            self._active_clusters.add(cname)
            self._conditions.clear()

            try:
                self._apply_clusters_to_adata()
                self._status_div.text = (
                    f"<p style='color:green;font-size:11px'>"
                    f"Cluster '{cname}' created.</p>"
                )
            except Exception as e:
                import traceback
                traceback.print_exc()
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            def _update():
                _render_conditions()
                self._rebuild_cluster_widgets(doc)
            doc.add_next_tick_callback(_update)

        def on_delete_cluster():
            cname = cluster_name_input.value.strip()
            if cname not in self._clusters:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>"
                    f"Cluster '{cname}' does not exist.</p>"
                )
                return

            del self._clusters[cname]
            self._active_clusters.discard(cname)

            if f"grp_{cname}" in self.sample.adata.obs.columns:
                self.sample.adata.obs.drop(columns=[f"grp_{cname}"], inplace=True)

            try:
                if self._clusters:
                    self._apply_clusters_to_adata()
                else:
                    # No clusters left — remove column from obs
                    if GROUPS_COL in self.sample.adata.obs.columns:
                        self.sample.adata.obs.drop(columns=[GROUPS_COL], inplace=True)

                self._status_div.text = (
                    f"<p style='color:#f0a500;font-size:11px'>"
                    f"Cluster '{cname}' deleted.</p>"
                )
            except Exception as e:
                import traceback
                traceback.print_exc()
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            doc.add_next_tick_callback(lambda: self._rebuild_cluster_widgets(doc))

        add_condition_btn.on_click(on_add_condition)
        clear_conditions_btn.on_click(on_clear_conditions)
        create_cluster_btn.on_click(on_create_cluster)
        delete_cluster_btn.on_click(on_delete_cluster)

        # ------------------------------------------------------------------
        # Layout
        # ------------------------------------------------------------------

        condition_row = row(gene_input, operator_select, value_input, logic_select, spacing=5)
        condition_btns = row(add_condition_btn, clear_conditions_btn, spacing=5)
        action_btns = row(create_cluster_btn, delete_cluster_btn, spacing=10)

        # --- DE section ---
        self._de_div = Div(text="", width=380)
        self._de_cluster_select = Select(
            title="Run DE for cluster",
            options=list(self._clusters.keys()) or ["—"],
            value=list(self._clusters.keys())[0] if self._clusters else "—",
            width=200
        )
        run_de_btn = Button(label="▶ Run DE", button_type="primary", width=120)

        def on_run_de():
            active = list(self._active_clusters)
            if not active:
                return
            if hasattr(self.component, '_figure_de_div'):
                self.component._figure_de_div.text = (
                    "<p style='color:#aaa;font-size:11px'>Running DE...</p>"
                )
            import threading
            def _run_in_thread():
                adata_copy = self.sample.adata.copy()
                sc.pp.log1p(adata_copy)
                results = {}
                for cname in active:
                    col_key = f"grp_{cname}"
                    if col_key not in adata_copy.obs.columns:
                        continue
                    in_cluster = adata_copy.obs[col_key].astype(bool)
                    n_in = in_cluster.sum()
                    n_out = (~in_cluster).sum()
                    if n_out == 0:
                        print(f"[DE] {cname}: todos los spots en cluster, saltando")
                        continue
                    if n_in < 3:
                        print(f"[DE] {cname}: muy pocos spots ({n_in}), saltando")
                        continue
                    adata_copy.obs['_de_group'] = adata_copy.obs[col_key].astype(bool).map(
                        {True: cname, False: 'Rest'}
                    ).astype(str)
                    try:
                        sc.tl.rank_genes_groups(
                            adata_copy, groupby='_de_group',
                            groups=[cname], reference='Rest',
                            method='wilcoxon', n_genes=adata_copy.n_vars,
                        )
                        results[cname] = sc.get.rank_genes_groups_df(adata_copy, group=cname)
                    except Exception as e:
                        print(f"[DE] Error en {cname}: {e}")
                        continue
                if not results:
                    return
                table_html = self._build_multi_de_table(results)
                plots_html = self._build_multi_plots_html(active, adata_copy, results)
                full_html = f"""
                <div style='display:flex;gap:15px;width:100%;height:600px'>
                    <div style='flex:0 0 400px;overflow-y:auto'>{table_html}</div>
                    <div style='flex:1;overflow-y:auto'>{plots_html}</div>
                </div>"""
                def _update():
                    if hasattr(self.component, '_figure_de_div'):
                        self.component._figure_de_div.text = full_html
                    if hasattr(self.component, '_metadata_panel'):
                        self.component._metadata_panel.update_de_table(
                            self.component, ', '.join(active), full_html
                        )
                doc.add_next_tick_callback(_update)

            thread = threading.Thread(target=_run_in_thread)
            thread.daemon = True
            thread.start()

        run_de_btn.on_click(on_run_de)

        return [
            Div(text="<b style='color:#888888'>Build conditions</b>"),
            condition_row,
            condition_btns,
            self._conditions_div,
            Div(text="<hr style='border-color:#444'>"),
            cluster_name_input,
            action_btns,
            self._status_div,
            Div(text="<b style='color:#888888'>Created clusters</b>"),
            self._cluster_layout,
            Div(text="<hr style='border-color:#444'>"),
            run_de_btn,
        ]

    # ------------------------------------------------------------------
    # Cluster widgets (checkbox list) — same pattern as ClusteringManager
    # ------------------------------------------------------------------

    def _build_cluster_widgets(self, doc) -> list:
        if not self._clusters:
            return [Div(
                text="<p style='color:#888;font-size:11px'>No clusters created yet.</p>"
            )]

        labels = []
        for cname, conditions in self._clusters.items():
            n_spots = self._count_spots_bool(cname)
            cond_str = " & ".join(
                f"{col}{op}{val}" for col, op, val, logic in conditions
            )
            labels.append(f"{cname}  ({n_spots} spots)  [{cond_str}]")

        self._checkbox_group = CheckboxGroup(
            labels=labels,
            active=[i for i, cname in enumerate(self._clusters)
                    if cname in self._active_clusters],
        )
        self._apply_visibility_mask()

        def _on_checkbox_change(attr, old, new):
            cluster_names = list(self._clusters.keys())
            self._active_clusters = {cluster_names[i] for i in new}
            doc.add_next_tick_callback(self._apply_visibility_mask)

        self._checkbox_group.on_change("active", _on_checkbox_change)
        return [self._checkbox_group]

    def _rebuild_cluster_widgets(self, doc):
        if self._cluster_layout is None:
            return
        self._cluster_layout.children = self._build_cluster_widgets(doc)
        # Actualizar opciones del selector DE
        if hasattr(self, '_de_cluster_select') and self._clusters:
            self._de_cluster_select.options = list(self._clusters.keys())
            self._de_cluster_select.value = list(self._clusters.keys())[0]


    def _count_spots_bool(self, cname: str) -> int:
        col_key = f"grp_{cname}"
        obs = self.sample.adata.obs
        if col_key not in obs.columns:
            return 0
        return int(obs[col_key].astype(bool).sum())

    # ------------------------------------------------------------------
    # Visibility mask — same pattern as ClusteringManager
    # ------------------------------------------------------------------

    def _apply_visibility_mask(self):
        if not self._clusters:
            return

        adata = self.sample.adata
        gdf = self.sample.gdf

        # 1. Calcular máscara de visibilidad
        visible = pd.Series([False] * len(gdf), index=gdf.index)
        for cname in self._active_clusters:
            col_key = f"grp_{cname}"
            if col_key in adata.obs.columns:
                cluster_series = self.sample.get_value(colname=col_key)
                if not isinstance(cluster_series, pd.Series):
                    cluster_series = pd.Series(cluster_series, index=gdf.index)
                cluster_series = cluster_series.reindex(gdf.index).fillna(False).astype(bool)
                visible |= cluster_series

        # 2. Reconstruir la columna custom_clusters
        col = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)
        for cname in self._clusters:
            if cname in self._active_clusters:
                col_key = f"grp_{cname}"
                if col_key in adata.obs.columns:
                    col[adata.obs[col_key].astype(bool)] = cname
        adata.obs[GROUPS_COL] = col

        # 3. Sincronizar la vista con el componente visual de Bokeh
        self.component._view = visible.tolist()
        
        # Si 'change_glyph_view' usa la API vieja, asegúrate de que haga esto internamente:
        # de lo contrario, fuerza la actualización aquí:
        if hasattr(self.component, 'source_view') and self.component.source_view is not None:
            from bokeh.models import IndexFilter
            # Convertimos la lista booleana a una lista de índices enteros (los que son True)
            indices_a_pintar = [i for i, x in enumerate(self.component._view) if x]
            self.component.source_view.filter = IndexFilter(indices=indices_a_pintar)

        # Invocar el refresco nativo del componente
        self.component.change_glyph_view()

    # ------------------------------------------------------------------
    # Apply clusters to adata.obs[GROUPS_COL]
    # ------------------------------------------------------------------

    def _apply_clusters_to_adata(self):
        from bokeh.palettes import Category20, Category10
        PALETTE = [*Category20[20], *Category10[10]]

        adata = self.sample.adata

        for cname, conditions in self._clusters.items():
            mask = self._build_mask(conditions, adata)
            adata.obs[f"grp_{cname}"] = mask.values
            if cname not in self._color_map:
                idx = len(self._color_map) % len(PALETTE)
                self._color_map[cname] = PALETTE[idx]
            print(f"[GroupsManager] {cname}: {mask.sum()} spots")

        self._rebuild_custom_clusters_col()

        # Crear/reemplazar parámetro color
        self.component.clustering_manager.apply_colors_to_column(GROUPS_COL)

        # Sobrescribir category_colors con nuestro _color_map fijo
        color_param = self.component.manager.get_param('color')
        if color_param is not None and hasattr(color_param, 'category_colors'):
            color_param.category_colors = dict(self._color_map)

            # Recalcular source con los colores correctos
            if GROUPS_COL in self.component.manager.cache.columns:
                self.component.manager.cache.drop(columns=[GROUPS_COL], inplace=True)

            obs_series = self.sample.adata.obs[GROUPS_COL]
            obs_names = pd.Series(self.component.manager.source.data['obs_name'])
            aligned = obs_names.map(obs_series)
            self.component.manager.cache[GROUPS_COL] = aligned.values

            print(f"[DEBUG custom_clusters] unique en adata.obs: {self.sample.adata.obs[GROUPS_COL].dropna().unique()}")
            print(f"[DEBUG aligned] unique: {pd.Series(aligned).dropna().unique()}")

            new_colors = color_param.get_values(
                pd.Series(self.component.manager.cache[GROUPS_COL]).copy()
            )
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            data['custom_clusters'] = list(aligned.fillna('—').values)
            self.component.manager.source.data = data

        # VALIDACIÓN — quitar después
        for cname in self._clusters:
            col_key = f"grp_{cname}"
            n_adata = int(self.sample.adata.obs[col_key].astype(bool).sum())
            n_custom = int((self.sample.adata.obs[GROUPS_COL] == cname).sum())
            print(f"[VALIDACIÓN] {cname}: grp_col={n_adata}, custom_clusters={n_custom}")

    def _rebuild_custom_clusters_col(self):
        """Recalcula custom_clusters — el último cluster activo gana en solapamientos."""
        adata = self.sample.adata
        col = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)

        for cname in self._clusters:  # orden de creación
            col_key = f"grp_{cname}"
            if col_key in adata.obs.columns:
                mask = adata.obs[col_key].astype(bool)
                col[mask] = cname  # último sobreescribe

        adata.obs[GROUPS_COL] = col
    # ------------------------------------------------------------------
    # Build boolean mask from conditions
    # ------------------------------------------------------------------

    def _build_mask(self, conditions, adata):
        obs = adata.obs
        mask = None

        for i, (col, op, val, logic) in enumerate(conditions):
            val_f = float(val)

            if col in adata.var_names:
                col_idx = list(adata.var_names).index(col)
                expr = adata.X[:, col_idx]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                series = pd.Series(expr, index=obs.index)
            elif col in obs.columns:
                series = obs[col].astype(float)
            else:
                print(f"[GroupsManager] Warning: '{col}' not found, skipping.")
                continue

            if op == ">":
                cond_mask = series > val_f
            elif op == "<":
                cond_mask = series < val_f
            elif op == ">=":
                cond_mask = series >= val_f
            elif op == "<=":
                cond_mask = series <= val_f
            else:
                cond_mask = series == val_f

            if mask is None:
                mask = cond_mask
            elif logic == "AND":
                mask = mask & cond_mask
            else:
                mask = mask | cond_mask

        if mask is None:
            return pd.Series([False] * len(obs), index=obs.index)

        return mask

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def _count_spots(self, cname: str) -> int:
        obs = self.sample.adata.obs
        if GROUPS_COL not in obs.columns:
            return 0
        return int((obs[GROUPS_COL] == cname).sum())
    
    def _run_de(self, cname: str) -> str:
        adata = self.sample.adata.copy()
        col_key = f"grp_{cname}"
        if col_key not in adata.obs.columns:
            return "<p style='color:red;padding:10px'>Cluster not found.</p>"
        
        adata.obs['_de_group'] = adata.obs[col_key].astype(bool).map(
            {True: cname, False: 'Rest'}
        ).astype(str)
        
        try:
            sc.tl.rank_genes_groups(
                adata, groupby='_de_group',
                groups=[cname], reference='Rest',
                method='wilcoxon', n_genes=adata.n_vars,
            )
        except Exception as e:
            return f"<p style='color:red;padding:10px'>DE error: {e}</p>"
        try: 
            result = sc.get.rank_genes_groups_df(adata, group=cname)
            result = result[['names', 'scores', 'pvals_adj', 'logfoldchanges']]
        except Exception as e:
            print(f"[plots] No hay DE para {cname}: {e}")
            return f"<p style='color:red;padding:10px'>No DE results for {cname}: {e}</p>"
        
        rows_html = ""
        for _, r in result.iterrows():
            lfc = float(r['logfoldchanges'])
            lfc_color = "#4ef1a1" if lfc >= 0 else "#f14e4e"
            rows_html += f"""<tr style='border-bottom:1px solid #333'>
                <td style='padding:5px 10px;color:#fff'>{r['names']}</td>
                <td style='padding:5px 10px;text-align:right;color:{lfc_color}'>{lfc:.2f}</td>
                <td style='padding:5px 10px;text-align:right;color:#aaa'>{float(r['pvals_adj']):.2e}</td>
                <td style='padding:5px 10px;text-align:right;color:#4e9af1'>{float(r['scores']):.2f}</td>
            </tr>"""
        
        return f"""
        <div style='padding:10px;font-size:13px;color:#ccc;width:100%;
                    height:400px;overflow-y:auto;'>
            <h4 style='color:#fff;margin:0 0 10px 0;position:sticky;top:0;
                    background:#1e1e2e;padding:5px 0'>DE: {cname} vs Rest</h4>
            <table style='border-collapse:collapse;font-size:13px;width:100%'>
                <thead style='position:sticky;top:32px;background:#1e1e2e'>
                    <tr style='border-bottom:1px solid #555'>
                        <th style='text-align:left;padding:5px 10px;color:#fff;font-weight:normal'>Gene</th>
                        <th style='text-align:right;padding:5px 10px;color:#fff;font-weight:normal'>L2FC</th>
                        <th style='text-align:right;padding:5px 10px;color:#fff;font-weight:normal'>Adj. p-val</th>
                        <th style='text-align:right;padding:5px 10px;color:#fff;font-weight:normal'>Score</th>
                    </tr>
                </thead>
                <tbody>{rows_html}</tbody>
            </table>
        </div>"""

    # ------------------------------------------------------------------
    # Implementación de los gráficos
    # ------------------------------------------------------------------

    def _build_plots_html(self, cname: str, adata,
                      cluster_color='#4ef1a1',
                      rest_color='#cccccc',
                      spatial_figsize=(4, 4),
                      violin_figsize=(3, 3.5),
                      dpi=100) -> str:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import io, base64

        col_key = f"grp_{cname}"
        obs = adata.obs.copy()
        in_cluster = obs[col_key].astype(bool)

        def fig_to_base64(fig):
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=100, bbox_inches='tight',
                        facecolor='#444')
            buf.seek(0)
            encoded = base64.b64encode(buf.read()).decode('utf-8')
            plt.close(fig)
            return encoded

        # --- Plot 1: Spatial scatter ---
        fig1, ax1 = plt.subplots(figsize=(4, 4))
        ax1.set_facecolor('#444')
        ax1.scatter(
            obs.loc[~in_cluster, 'pxl_col_in_fullres'],
            -obs.loc[~in_cluster, 'pxl_row_in_fullres'],
            c='#cccccc', s=2, alpha=0.3, linewidths=0
        )
        ax1.scatter(
            obs.loc[in_cluster, 'pxl_col_in_fullres'],
            -obs.loc[in_cluster, 'pxl_row_in_fullres'],
            c='#4ef1a1', s=4, alpha=1.0, linewidths=0
        )
        ax1.set_title(f'Spatial: {cname}', color='white', fontsize=11)
        ax1.axis('off')
        img1 = fig_to_base64(fig1)

        # --- Plot 2: Strip plot top genes ---
        import scanpy as sc
        result = sc.get.rank_genes_groups_df(adata, group=cname)
        top_genes = result[result['pvals_adj'] < 0.05].head(5)['names'].tolist()
        if not top_genes:
            top_genes = result.head(5)['names'].tolist()

        n = len(top_genes)
        fig2, axes = plt.subplots(1, n, figsize=(3 * n, 3.5))
        if n == 1:
            axes = [axes]
        fig2.patch.set_facecolor('#444')

        for ax, gene in zip(axes, top_genes):
            ax.set_facecolor('#444')
            if gene not in adata.var_names:
                continue
            idx = list(adata.var_names).index(gene)
            expr = adata.X[:, idx]
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()

            in_vals = expr[in_cluster.values]
            out_vals = expr[~in_cluster.values]

            # Strip plot con jitter
            import numpy as np
            jitter = 0.15
            ax.scatter(
                np.random.uniform(1 - jitter, 1 + jitter, len(in_vals)),
                in_vals, c='#4ef1a1', s=3, alpha=0.5, linewidths=0
            )
            ax.scatter(
                np.random.uniform(2 - jitter, 2 + jitter, len(out_vals)),
                out_vals, c='#aaaaaa', s=3, alpha=0.3, linewidths=0
            )
            ax.set_title(gene, color='white', fontsize=10)
            ax.set_xticks([1, 2])
            ax.set_xticklabels([cname, 'Rest'], color='#aaaaaa', fontsize=8)
            ax.tick_params(colors='#aaaaaa')
            ax.spines['bottom'].set_color('#444')
            ax.spines['left'].set_color('#444')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            for label in ax.get_yticklabels():
                label.set_color('#aaaaaa')

        plt.tight_layout()
        img2 = fig_to_base64(fig2)

        return f"""
        <img src='data:image/png;base64,{img1}' style='width:100%;margin-bottom:10px'/>
        <img src='data:image/png;base64,{img2}' style='width:100%;margin-bottom:10px'/>
        """
    def _build_multi_de_table(self, results: dict) -> str:
        # Unir todos los genes
        import pandas as pd
        all_genes = set()
        for df in results.values():
            all_genes.update(df['names'].tolist())
        
        # Construir tabla unificada
        merged = pd.DataFrame({'Gene': sorted(all_genes)})
        for cname, df in results.items():
            df = df.set_index('names')
            color = self._color_map.get(cname, '#4ef1a1')
            merged[f'{cname}_lfc'] = merged['Gene'].map(df['logfoldchanges']).round(2)
            merged[f'{cname}_pval'] = merged['Gene'].map(df['pvals_adj'])
        
        # Ordenar por score del primer cluster
        first = list(results.keys())[0]
        first_df = results[first].set_index('names')
        merged['_sort'] = merged['Gene'].map(first_df['scores']).fillna(0)
        merged = merged.sort_values('_sort', ascending=False).drop(columns='_sort')

        # Headers
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
                lfc = row.get(f'{cname}_lfc', float('nan'))
                pval = row.get(f'{cname}_pval', float('nan'))
                try:
                    lfc_color = '#4ef1a1' if float(lfc) >= 0 else '#f14e4e'
                    lfc_str = f"{float(lfc):.2f}"
                except:
                    lfc_color = '#888'
                    lfc_str = '—'
                try:
                    pval_str = f"{float(pval):.2e}"
                    pval_color = '#fff' if float(pval) < 0.05 else '#666'
                except:
                    pval_str = '—'
                    pval_color = '#666'
                cells += f"""
                    <td style='padding:4px 8px;text-align:right;color:{lfc_color}'>{lfc_str}</td>
                    <td style='padding:4px 8px;text-align:right;color:{pval_color}'>{pval_str}</td>"""
            rows_html += f"<tr style='border-bottom:1px solid #222'>{cells}</tr>"

        return f"""
        <div style='font-size:12px;height:100%'>
            <h4 style='color:#fff;margin:0 0 8px 0;position:sticky;top:0;
                    background:#000;padding:5px 0'>DE vs Rest</h4>
            <table style='border-collapse:collapse;width:100%'>
                <thead style='position:sticky;top:32px;background:#000'>
                    <tr><th style='padding:5px 10px;color:#888;text-align:left'>Gene</th>
                        {cluster_headers}</tr>
                    <tr><th></th>{sub_headers}</tr>
                </thead>
                <tbody>{rows_html}</tbody>
            </table>
        </div>"""

    def _build_multi_plots_html(self, active_clusters: list, adata, results: dict) -> str:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import io, base64, numpy as np

        def fig_to_base64(fig):
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=100, bbox_inches='tight', facecolor='#000000')
            buf.seek(0)
            enc = base64.b64encode(buf.read()).decode('utf-8')
            plt.close(fig)
            return enc

        obs = adata.obs.copy()

        # --- Spatial: todos los clusters ---
        fig1, ax1 = plt.subplots(figsize=(4, 4))
        ax1.set_facecolor('#000000')
        # Fondo gris para todos
        ax1.scatter(
            obs['pxl_col_in_fullres'], -obs['pxl_row_in_fullres'],
            c='#333333', s=2, alpha=0.3, linewidths=0
        )
        for cname in active_clusters:
            col_key = f"grp_{cname}"
            if col_key not in obs.columns:
                continue
            mask = obs[col_key].astype(bool)
            color = self._color_map.get(cname, '#4ef1a1')
            ax1.scatter(
                obs.loc[mask, 'pxl_col_in_fullres'],
                -obs.loc[mask, 'pxl_row_in_fullres'],
                c=color, s=5, alpha=1.0, linewidths=0, label=cname
            )
        ax1.legend(facecolor='#111', labelcolor='white', fontsize=8,
                loc='lower right', framealpha=0.7)
        ax1.set_title('Spatial', color='white', fontsize=11)
        ax1.axis('off')
        img_spatial = fig_to_base64(fig1)

        # --- Strip plots top genes por cluster ---
        strips_html = ""
        for cname in active_clusters:
            col_key = f"grp_{cname}"
            if col_key not in obs.columns:
                continue
            in_cluster = obs[col_key].astype(bool)
            color = self._color_map.get(cname, '#4ef1a1')

            # Top 5 genes significativos
            import scanpy as sc
            if cname not in results:
                continue
            result = results[cname]
            top_genes = result[result['pvals_adj'] < 0.05].head(5)['names'].tolist()
            if not top_genes:
                top_genes = result.head(5)['names'].tolist()

            n = len(top_genes)
            if n == 0:
                continue

            fig2, axes = plt.subplots(1, n, figsize=(2.5 * n, 3))
            if n == 1:
                axes = [axes]
            fig2.patch.set_facecolor('#000000')

            for ax, gene in zip(axes, top_genes):
                ax.set_facecolor('#000000')
                if gene not in adata.var_names:
                    continue
                idx = list(adata.var_names).index(gene)
                expr = adata.X[:, idx]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = expr.flatten()

                jitter = 0.15
                ax.scatter(
                    np.random.uniform(1-jitter, 1+jitter, in_cluster.sum()),
                    expr[in_cluster.values],
                    c=color, s=3, alpha=0.6, linewidths=0
                )
                ax.scatter(
                    np.random.uniform(2-jitter, 2+jitter, (~in_cluster).sum()),
                    expr[~in_cluster.values],
                    c='#555555', s=2, alpha=0.3, linewidths=0
                )
                ax.set_title(gene, color='white', fontsize=9)
                ax.set_xticks([1, 2])
                ax.set_xticklabels([cname[:8], 'Rest'], color='#aaaaaa', fontsize=7)
                ax.tick_params(colors='#aaaaaa')
                for spine in ['bottom', 'left']:
                    ax.spines[spine].set_color('#444')
                for spine in ['top', 'right']:
                    ax.spines[spine].set_visible(False)
                for label in ax.get_yticklabels():
                    label.set_color('#aaaaaa')

            plt.tight_layout()
            img_strip = fig_to_base64(fig2)
            strips_html += f"""
            <p style='color:#aaa;font-size:11px;margin:8px 0 4px 0'>
                <span style='color:{color}'>●</span> {cname} — top genes
            </p>
            <img src='data:image/png;base64,{img_strip}' style='width:100%'/>"""

        return f"""
        <img src='data:image/png;base64,{img_spatial}' style='width:100%;margin-bottom:12px'/>
        {strips_html}"""
