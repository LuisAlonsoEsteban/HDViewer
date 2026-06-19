import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.models import (
    Button, Div, Select, TextInput,
    AutocompleteInput, CheckboxGroup,
)
from bokeh.layouts import column, row

from bokeh.models import Button, Div, Select
from bokeh.models import RadioButtonGroup
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
        self._de_merged = None  # guardamos el DataFrame para el buscador
        self._de_results = {}   # guardamos results para regenerar
        self._de_table_div = None  # referencia al Div de la tabla
        self._last_plots_html = ""

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

        de_mode = RadioButtonGroup(labels=["vs Rest", "Entre clusters"], active=0)

        self._status_div = Div(text="", width=380)

        # Cluster list with checkboxes — updated in-place like ClusteringManager
        if self._cluster_layout is None:
            self._cluster_layout = column(
                name="groups_cluster_controls",
                children=self._build_cluster_widgets(doc),
            )
        # --- DE section ---
        from bokeh.models import TextInput as BokehTextInput

        self._de_table_div = Div(text="", width=380, height=500)
        de_search_input = BokehTextInput(
            placeholder="Search gene...",
            width=380
        )

        def on_de_search(attr, old, new):
            if self._de_merged is None or not self._de_results:
                return
            
            if new:
                filtered = self._de_merged[
                    self._de_merged['Gene'].str.lower().str.contains(new.lower())
                ]
            else:
                filtered = self._de_merged.head(200)
            
            def _update():
                table_html = self._render_de_table_html(filtered, self._de_results, len(self._de_merged))
                plots_html = getattr(self, '_last_plots_html', '')
                
                full_html = f"""
                <div style='display:grid;grid-template-columns:380px 1fr;gap:15px;width:100%'>
                    <div style='overflow-y:auto;max-height:600px'>{table_html}</div>
                    <div style='overflow-y:auto;max-height:600px'>{plots_html}</div>
                </div>"""
                
                if hasattr(self.component, '_figure_de_div'):
                    self.component._figure_de_div.text = full_html
            doc.add_next_tick_callback(_update)

        de_search_input.on_change('value', on_de_search)

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

            mode = de_mode.active
            if hasattr(self.component, '_figure_de_div'):
                self.component._figure_de_div.text = (
                    "<p style='color:#aaa;font-size:11px'>Running DE...</p>"
                )
            import threading
            def _run_in_thread():
                if getattr(self, '_de_adata', None) is None:
                    _a = self.sample.adata.copy()
                    sc.pp.filter_genes(_a, min_cells=3)
                    sc.pp.log1p(_a)  # elimino genes casi vacíos
                    self._de_adata = _a
                adata_copy = self._de_adata
                # refrescar columnas de clusters (pueden haberse creado/borrado)
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
                                method='wilcoxon', 
                            )
                            results[cname] = sc.get.rank_genes_groups_df(adata_copy, group=cname)
                        except Exception as e:
                            print(f"[DE] Error en {cname}: {e}")
                            continue
                else:
                    # --- Entre clusters seleccionados ---
                    obs_grp = adata_copy.obs[GROUPS_COL].astype(str)
                    in_sel = adata_copy.obs[GROUPS_COL].notna() & obs_grp.isin(active)
                    sub = adata_copy[in_sel.values].copy()
                    sub.obs['_de_group'] = sub.obs[GROUPS_COL].astype(str)
                    plot_adata = sub

                    counts = sub.obs['_de_group'].value_counts()
                    usable = [c for c in active if counts.get(c, 0) >= 3]
                    if len(usable) < 2:
                        print("[DE] Se necesitan al menos 2 clusters con ≥3 bins para comparar entre ellos")
                        return

                    sc.tl.rank_genes_groups(
                        sub, groupby='_de_group',
                        groups=usable, reference='rest',
                        method='wilcoxon',
                    )
                    for cname in usable:
                        results[cname] = sc.get.rank_genes_groups_df(sub, group=cname)
                if not results:
                    return

                if mode != 0:
                    _keys = list(results.keys())
                    self._de_header = (f"DE: {_keys[0]} vs {_keys[1]}"
                                       if len(_keys) == 2
                                       else "DE entre clusters seleccionados")
                else:
                    self._de_header = "DE vs Rest"
                table_html = self._build_multi_de_table(results)
                plots_html = self._build_multi_plots_html(active, plot_adata, results, mode)
                self._last_plots_html = plots_html

                print("table len =", len(table_html))
                print("plots len =", len(plots_html))

                full_html = f"""
                <div style='display:grid;grid-template-columns:380px 1fr;gap:15px;width:100%'>
                    <div style='overflow-y:auto;max-height:600px'>{table_html}</div>
                    <div style='overflow-y:auto;max-height:600px'>{plots_html}</div>
                </div>"""

                print(f"[_run_in_thread] full_html len={len(full_html)}")

                def _update():
                    print("HAS IMG =", "<img" in full_html)
                    print("HAS BASE64 =", "base64" in full_html)
                    print(f"[_update] ejecutándose, html len={len(full_html)}")
                    if hasattr(self.component, '_figure_de_div'):
                        self.component._figure_de_div.text = f"""
                        <div style='display:grid;grid-template-columns:380px 1fr;gap:15px;width:100%'>
                            <div style='overflow-y:auto;max-height:600px'>{table_html}</div>
                            <div style='overflow-y:auto;max-height:600px'>{plots_html}</div>
                        </div>"""

                doc.add_next_tick_callback(_update)

            thread = threading.Thread(target=_run_in_thread)
            thread.daemon = True
            thread.start()

        run_de_btn.on_click(on_run_de)

        _alias = getattr(self.component, "alias", None)
        _geom_name = (str(_alias) + "s") if _alias else "N/A"
        _geom_label_div = Div(text="<b style='color:#888888'>Viewing:</b> " + _geom_name)
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
            de_mode,
            run_de_btn,
            _geom_label_div,
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
        try:
            self.component.clustering_manager.apply_colors_to_column(GROUPS_COL)
        except ValueError:
            color_param = self.component.manager.get_param('color')
            if color_param is not None:
                color_param.colname = GROUPS_COL
                self.component.manager.update_param('color')
        # Sobrescribir category_colors con nuestro _color_map fijo
        color_param = self.component.manager.get_param('color')
        if color_param is not None and hasattr(color_param, 'category_colors'):
            color_param.category_colors = dict(self._color_map)
            if hasattr(self.component, 'loader'):
                self.component.manager.update_param('alpha', avoid_scaling=True)
                self.component.loader.force_reload()
                return
            obs_series = self.sample.adata.obs[GROUPS_COL]
            obs_names = pd.Series(self.component.manager.source.data['obs_name'])
            aligned = obs_names.map(obs_series)
            #if GROUPS_COL not in self.component.manager.cache.columns:
            self.component.manager.cache[GROUPS_COL] = aligned.values
            new_colors = color_param.get_values(
                pd.Series(self.component.manager.cache[GROUPS_COL]).copy()
            )
            data = dict(self.component.manager.source.data)
            data['color'] = new_colors
            data['custom_clusters'] = ['—' if pd.isna(c) else c for c in aligned]
            self.component.manager.source.data = data
            self.component.manager.fetch_data()
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

        adata.obs[GROUPS_COL] = col.astype('category')
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
            parts = ax.violinplot(
                [expr[in_cluster.values], expr[~in_cluster.values]],
                positions=[1, 2],
                showmedians=True,
                showextrema=False
            )
            # Color del cluster
            parts['bodies'][0].set_facecolor(color)
            parts['bodies'][0].set_alpha(0.7)
            # Color del resto
            parts['bodies'][1].set_facecolor('#555555')
            parts['bodies'][1].set_alpha(0.5)
            # Línea de mediana
            parts['cmedians'].set_color('white')
            parts['cmedians'].set_linewidth(1.5)

            ax.set_xticks([1, 2])
            ax.set_xticklabels([cname, 'Rest'], color='#aaaaaa', fontsize=7)

        plt.tight_layout()
        img2 = fig_to_base64(fig2)

        return f"""
        <img src='data:image/png;base64,{img1}' style='width:100%;margin-bottom:10px'/>
        <img src='data:image/png;base64,{img2}' style='width:100%;margin-bottom:10px'/>
        """
    def _build_multi_de_table(self, results: dict, filter_text: str = "") -> str:
        import pandas as pd
        all_genes = set()
        for df in results.values():
            all_genes.update(df['names'].tolist())

        merged = pd.DataFrame({'Gene': sorted(all_genes)})
        for cname, df in results.items():
            df = df.set_index('names')
            merged[f'{cname}_lfc'] = merged['Gene'].map(df['logfoldchanges']).round(2)
            merged[f'{cname}_pval'] = merged['Gene'].map(df['pvals_adj'])

        first = list(results.keys())[0]
        first_df = results[first].set_index('names')
        merged['_sort'] = merged['Gene'].map(first_df['scores']).fillna(0)
        merged = merged.sort_values('_sort', ascending=False).drop(columns='_sort')

        # Guardar para el buscador
        self._de_merged = merged
        self._de_results = results

        # Filtrar si hay texto
        if filter_text:
            merged = merged[merged['Gene'].str.lower().str.contains(filter_text.lower())]
        else:
            merged = merged.head(200)  # limitar a 50 por defecto

        return self._render_de_table_html(merged, results)

    def _render_de_table_html(self, merged, results, n_total=None) -> str:
        n_shown = len(merged)
        if n_total is None:
            n_total = n_shown
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
            rows_html += f"<tr class='de-row' data-gene='{row['Gene'].lower()}' style='border-bottom:1px solid #222'>{cells}</tr>"

        n_shown = len(merged)
        n_total = len(self._de_merged) if self._de_merged is not None else n_shown
        footer = f"<p style='color:#666;font-size:11px;margin:8px 0 0 0'>Showing {n_shown} of {n_total} genes</p>"

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
 
    def _build_multi_plots_html(self, active_clusters: list, adata, results: dict, mode: int = 0) -> str:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import io, base64, numpy as np
        import uuid, os

        sample_type = getattr(self.component, 'sample_type', 'spots')
        print("COMPONENT ALIAS:", getattr(self.component, 'alias', 'UNKNOWN'))
        print("SAMPLE TYPE:", sample_type)
        if sample_type is None:
            raise ValueError(f"Component has no sample_type: {self.component}")

        def fig_to_url(fig):
            import pathlib
            import uuid
            filename = f"de_plot_{uuid.uuid4().hex}.png"
            filepath = pathlib.Path("/scratch/tcondegerma/VisiumVisualizer/bokeh_server/static/data") / filename
            fig.savefig(str(filepath), format='png', dpi=100,
                        bbox_inches='tight', facecolor='#000000')
            plt.close(fig)
            return f"/bokeh_server/static/data/{filename}?t={uuid.uuid4().hex[:6]}"

        obs = adata.obs.copy()
        n_genes = 5

        # Recopilar top genes por cluster
        cluster_genes = {}
        for cname in active_clusters:
            if cname not in results:
                continue
            result = results[cname]
            sig = result[(result['pvals_adj'] < 0.05) & (result['logfoldchanges'] > 0)].sort_values('scores', ascending=False)
            top = sig.head(n_genes)['names'].tolist()
            if len(top) < n_genes:
                rest = result[(~result['names'].isin(top)) & (result['logfoldchanges'] > 0)].sort_values('scores', ascending=False)
                top += rest.head(n_genes - len(top))['names'].tolist()
            cluster_genes[cname] = top
        violins_html = ""

        print("ACTIVE CLUSTERS:", active_clusters)
        print("RESULTS KEYS:", list(results.keys()))
        print("SAMPLE TYPE:", sample_type)

        # Paleta de colores distintivos para los genes
        GENE_COLORS = [
            '#4ef1a1', '#f1a14e', '#4e9af1', '#f14e4e', '#c44ef1',
            '#f1e44e', '#4ef1e4', '#f14ebc', '#8af14e', '#f17a4e'
        ]

        for cname in active_clusters:
            if cname not in cluster_genes or not cluster_genes[cname]:
                continue
            if mode != 0 and '_de_group' in obs.columns:
                in_cluster = (obs['_de_group'].astype(str) == cname)
                _others = [c for c in active_clusters if c != cname]
                rest_label = " + ".join(_others) if _others else "Otros"
            else:
                col_key = f"grp_{cname}"
                if col_key not in obs.columns:
                    continue
                in_cluster = obs[col_key].astype(bool)
                rest_label = "Rest"
            print(f"[VIOLIN DEBUG] cname={cname} mode={mode} has_de_group={'_de_group' in obs.columns} rest_label={rest_label!r}")
            cluster_color = self._color_map.get(cname, '#4ef1a1')
            top_genes = cluster_genes[cname]
            n = len(top_genes)

            if sample_type == 'spots':
                dot_means_in, dot_pct_in = [], []
                dot_means_out, dot_pct_out = [], []
                valid_genes = []

                for gene in top_genes:
                    if gene not in adata.var_names:
                        continue
                    idx = list(adata.var_names).index(gene)
                    expr = adata.X[:, idx]
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    in_vals = np.log1p(expr[in_cluster.values])
                    out_vals = np.log1p(expr[~in_cluster.values])
                    dot_means_in.append(np.nanmean(in_vals))
                    dot_means_out.append(np.nanmean(out_vals))
                    dot_pct_in.append(np.mean(in_vals > 0) * 100)
                    dot_pct_out.append(np.mean(out_vals > 0) * 100)
                    valid_genes.append(gene)

                if not valid_genes:
                    continue

                lfc_map = dict(zip(results[cname]['names'], results[cname]['logfoldchanges']))
                lfc_values = [lfc_map.get(g, 0) for g in valid_genes]
                order = np.argsort([lfc_map.get(g, 0) for g in valid_genes])[::-1]
                dot_means_in  = [lfc_values[i]  for i in order]
                dot_means_out = [lfc_values[i] for i in order]
                dot_pct_in    = [dot_pct_in[i]    for i in order]
                dot_pct_out   = [dot_pct_out[i]   for i in order]
                valid_genes   = [valid_genes[i]   for i in order]
                lfc_values    = [lfc_values[i]    for i in order]
                n_valid = len(valid_genes)
                gene_colors = [GENE_COLORS[i % len(GENE_COLORS)] for i in range(n_valid)]

                fig_h = max(3.0, 0.3 * n_valid + 2.0)
                fig, ax = plt.subplots(figsize=(max(7, 0.9 * n_valid + 3), 5))
                fig.patch.set_facecolor('#000000')
                ax.set_facecolor('#0a0a0a')

                y_pos = np.arange(n_valid)
                max_size = 400

                # Banda de fondo alternada por gen — ahora en X
                for i, x in enumerate(y_pos):
                    bg = '#141414' if i % 2 == 0 else '#0a0a0a'
                    ax.axvspan(x - 0.5, x + 0.5, color=bg, zorder=0)

                # Línea horizontal en y=0
                ax.axhline(0, color='#333333', linewidth=0.8, zorder=1)

                # Puntos
                for i, (gene, gc) in enumerate(zip(valid_genes, gene_colors)):
                    x = y_pos[i]
                    s_in = max_size * (dot_pct_in[i] / 100)
                    s_out = max_size * (dot_pct_out[i] / 100)

                    if lfc_values[i] > 0 or dot_pct_out[i] > 0:
                        ax.scatter(
                            x - 0.18, lfc_values[i],
                            s=max(s_out, 15), c='#444444',
                            alpha=0.7, linewidths=0.8,
                            edgecolors='#888888', zorder=3
                        )
                        ax.plot(
                            [x - 0.18, x + 0.18],
                            [lfc_values[i], lfc_values[i]],
                            color='#333333', linewidth=0.6, zorder=2
                        )
                    # Punto cluster
                    ax.scatter(
                        x + 0.18, lfc_values[i],
                        s=max(s_in, 15), c=gc,
                        alpha=0.9, linewidths=0.8,
                        edgecolors='white', zorder=4
                    )

                ax.set_xticks(y_pos)
                ax.set_xticklabels(valid_genes, fontsize=12, rotation=30, ha='right')
                for tick, gc in zip(ax.get_xticklabels(), gene_colors):
                    tick.set_color(gc)

                ax.set_ylabel('Log2 Fold Change', color='white', fontsize=11)
                ax.set_xlabel('')
                ax.tick_params(colors='white', labelsize=7)
                ax.tick_params(axis='x', labelsize=7)
                plt.setp(ax.get_xticklabels(), rotation=20, ha='right')
                ax.set_ylim(bottom=-0.05)
                ax.set_xlim(-0.6, n_valid - 0.4)

                for spine in ['bottom', 'left']:
                    ax.spines[spine].set_color('#333')
                for spine in ['top', 'right']:
                    ax.spines[spine].set_visible(False)
                for label in ax.get_yticklabels():
                    label.set_color('white')

                # Leyenda grupo (cluster vs rest)
                from matplotlib.lines import Line2D
                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#aaaaaa',
                        markersize=8, label=f'{cname} (colored)', linestyle='None'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#444444',
                        markersize=6, label=rest_label, linestyle='None',
                        markeredgecolor='#888888'),
                ]
        
                leg1 = ax.legend(
                    handles=legend_elements,
                    loc='upper left', bbox_to_anchor=(0.9, 0.99),
                    fontsize=10, framealpha=0.3, labelcolor='white',
                    facecolor='#111111', edgecolor='#444444',
                    title='Group', title_fontsize=10,
                    borderaxespad=0
                )
                plt.setp(leg1.get_title(), color='#aaaaaa')
                ax.add_artist(leg1)

                # Leyenda tamaño de puntos
                pct_levels = [25, 50, 100]
                size_handles = [
                    plt.scatter([], [], s=max_size * p / 100,
                                c='#888888', alpha=0.7,
                                edgecolors='white', linewidths=0.5,
                                label=f'{p}%')
                    for p in pct_levels
                ]
                
                leg2 = ax.legend(
                    handles=size_handles,
                    title='% expressing', title_fontsize=10,
                    fontsize=10, loc='upper left', bbox_to_anchor=(1.02, 0.5),
                    framealpha=0.3, labelcolor='white',
                    facecolor='#111111', edgecolor='#444444',
                    handletextpad=0.3, borderpad=0.6,
                    borderaxespad=0
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
                </div>"""

            else:
                # --- Violin para bins y cells ---
                import pandas as pd
                import seaborn as sns

                fig_w = max(2.5 * n, 5)
                fig, axes = plt.subplots(1, n, figsize=(fig_w, 3.5))
                if n == 1:
                    axes = [axes]
                fig.patch.set_facecolor('#000000')

                for ax, gene in zip(axes, top_genes):
                    ax.set_facecolor('#111111')
                    if gene not in adata.var_names:
                        continue
                    idx = list(adata.var_names).index(gene)
                    expr = adata.X[:, idx]
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()
                    else:
                        expr = np.asarray(expr).flatten()
                    in_vals = np.log1p(expr[in_cluster.values])
                    out_vals = np.log1p(expr[~in_cluster.values])

                    in_nz = in_vals[in_vals > 0]
                    out_nz = out_vals[out_vals > 0]
                    if len(in_nz) < 3:
                        in_nz = in_vals
                    if len(out_nz) < 3:
                        out_nz = out_vals

                    df_v = pd.DataFrame({
                        'Expression': np.concatenate([in_nz, out_nz]),
                        'Group': [cname] * len(in_nz) + [rest_label] * len(out_nz)
                    })
                    sns.violinplot(
                        data=df_v, x='Group', y='Expression', ax=ax,
                        cut=0, inner='quartile',
                        hue='Group', legend=False,
                        palette={cname: cluster_color, rest_label: '#555555'},
                        order=[cname, rest_label],
                    )
                    print('DF GROUPS:', list(df_v['Group'].unique()))
                    #ymax = np.nanmax([np.nanmean(in_vals), np.nanmean(out_vals)])
                    #if np.isnan(ymax) or ymax <= 0:
                    #    ymax = 1
                    all_vals = np.concatenate([in_nz, out_nz])
                    ymax = np.nanpercentile(all_vals, 95) if len(all_vals) > 0 else 1.0
                    if np.isnan(ymax) or ymax <= 0:
                        ymax = 1.0
                    ax.set_ylim(0, ymax * 1.3)
                    ax.set_ylim(0, ymax * 1.4)
                    ax.set_title(gene, color='white', fontsize=9)
                    ax.set_xlabel('')
                    ax.tick_params(colors='#aaaaaa', labelsize=10)
                    for spine in ['bottom', 'left']:
                        ax.spines[spine].set_color('#444')
                    for spine in ['top', 'right']:
                        ax.spines[spine].set_visible(False)
                    for label in ax.get_xticklabels() + ax.get_yticklabels():
                        label.set_color('#aaaaaa')

                plt.tight_layout()
                img_url = fig_to_url(fig)
                img_w = max(250 * n, 500)
                violins_html += f"""
                <div style='margin-bottom:16px'>
                    <p style='color:#aaa;font-size:11px;margin:0 0 4px 0'>
                        <span style='color:{cluster_color}'>●</span>
                        <b style='color:{cluster_color}'>{cname}</b> — top {n} genes
                    </p>
                    <div style='overflow-x:auto;width:100%'>
                        <img src='{img_url}' style='width:{img_w}px;max-width:none;display:block'/>
                    </div>
                </div>"""

        if not violins_html:
            violins_html = "<p style='color:#888;font-size:12px'>No results to display.</p>"

        return violins_html