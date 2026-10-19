import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.models import (
    Button, Div, Select, TextInput,
    AutocompleteInput, CheckboxGroup,
    RadioButtonGroup,
)
from bokeh.layouts import column, row

from .groups_de import DEMixin

OPERATORS  = [">", "<", ">=", "<=", "=="]
GROUPS_COL = "custom_clusters"


class GroupsManager(DEMixin):
    """
    Manages creation of custom clusters based on gene/metadata conditions.
    Stores results in adata.obs[GROUPS_COL], same structure as leiden.

    Each spot gets the name of the first matching cluster, or None if unassigned.
    """

    def __init__(self, component):
        self.component = component
        self._clusters: dict        = {}   # cluster_name → [conditions]
        self._conditions: list      = []   # conditions being built [(col, op, val, logic)]
        self._active_clusters: set  = set()
        self._checkbox_group        = None
        self._cluster_layout        = None
        self._color_map: dict       = {}
        self._conditions_div        = None
        self._status_div            = None
        self._de_merged             = None
        self._de_results: dict      = {}
        # DE has its own selection: hiding a group on the map does not change it
        self._de_selected: set      = set()
        self._de_groups_checkbox    = None
        self._de_status_div         = None
        self._last_plots_html: str  = ""

    @property
    def sample(self):
        return self.component.manager.sample

    @property
    def is_clustered(self) -> bool:
        return (
            GROUPS_COL in self.sample.adata.obs.columns
            and len(self._clusters) > 0
        )

    # ------------------------------------------------------------------
    # Public: get_widgets
    # ------------------------------------------------------------------

    def get_widgets(self, doc) -> list:
        colnames = self.sample.get_colnames()

        gene_input = AutocompleteInput(
            title="Gene / Metadata",
            value=colnames[0] if colnames else "",
            completions=colnames,
            min_characters=1,
            case_sensitive=False,
            restrict=False,
            width=235,
        )
        operator_select = Select(title="Operator", options=OPERATORS, value=">", width=70)
        value_input     = TextInput(title="Value", value="1", width=80)
        logic_select    = Select(title="Logic", options=["AND", "OR"], value="AND", width=70)

        self._conditions_div = Div(
            text="<p style='color:#888;font-size:11px'>No conditions added yet.</p>",
            width=270,
        )
        cluster_name_input = TextInput(value="Cluster A", title="Cluster name", width=200)

        add_condition_btn    = Button(label="+ Add condition",  button_type="success", width=150)
        clear_conditions_btn = Button(label="Clear",            button_type="warning",  width=80)
        create_cluster_btn   = Button(label="✓ Create cluster", button_type="success", width=130)
        delete_cluster_btn   = Button(label="🗑 Delete cluster", button_type="danger",  width=130)

        self._status_div = Div(text="", width=270)

        if self._cluster_layout is None:
            self._cluster_layout = column(
                name="groups_cluster_controls",
                children=self._build_cluster_widgets(doc),
            )

        # --- Render helpers ---

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

        # --- Callbacks ---

        def on_add_condition():
            col   = gene_input.value.strip()
            op    = operator_select.value
            val   = value_input.value.strip()
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

            if col in self.sample.adata.var_names:
                var_idx = {g: i for i, g in enumerate(self.sample.adata.var_names)}
                expr = self.sample.adata.X[:, var_idx[col]]
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
                    f"<p style='color:red;font-size:11px'>Cluster '{cname}' already exists.</p>"
                )
                return

            self._clusters[cname] = list(self._conditions)
            self._active_clusters.add(cname)
            self._de_selected.add(cname)
            self._conditions.clear()

            try:
                self._apply_clusters_to_adata()
                self._status_div.text = (
                    f"<p style='color:green;font-size:11px'>Cluster '{cname}' created.</p>"
                )
            except Exception as e:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            def _update():
                _render_conditions()
                self._rebuild_cluster_widgets(doc)
                self._refresh_de_groups()

            doc.add_next_tick_callback(_update)

        def on_delete_cluster():
            cname = cluster_name_input.value.strip()
            if cname not in self._clusters:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Cluster '{cname}' does not exist.</p>"
                )
                return

            del self._clusters[cname]
            self._active_clusters.discard(cname)
            self._de_selected.discard(cname)
            self._clear_de_if_cluster_deleted()

            col_key = f"grp_{cname}"
            if col_key in self.sample.adata.obs.columns:
                self.sample.adata.obs.drop(columns=[col_key], inplace=True)

            try:
                if self._clusters:
                    self._apply_clusters_to_adata()
                else:
                    if GROUPS_COL in self.sample.adata.obs.columns:
                        self.sample.adata.obs.drop(columns=[GROUPS_COL], inplace=True)
                    if hasattr(self.component, 'loader'):
                        # Bins: quitar la máscara de grupos del view_filter.
                        mgr = self.component.manager
                        mgr.view_filter = np.ones(len(mgr.indexes), dtype=bool)
                        self.component._view = mgr.view_filter.tolist()
                        self.component.change_glyph_view()
                self._status_div.text = (
                    f"<p style='color:#f0a500;font-size:11px'>Cluster '{cname}' deleted.</p>"
                )
            except Exception as e:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            def _update():
                self._rebuild_cluster_widgets(doc)
                self._refresh_de_groups()
            doc.add_next_tick_callback(_update)

        add_condition_btn.on_click(on_add_condition)
        clear_conditions_btn.on_click(on_clear_conditions)
        create_cluster_btn.on_click(on_create_cluster)
        delete_cluster_btn.on_click(on_delete_cluster)

        # --- Layout ---
        # Gene en su propia fila y operador/valor/lógica debajo: todo en una
        # fila (~400 px) no cabía en el panel izquierdo y quedaba cortado.
        condition_row  = column(
            gene_input,
            row(operator_select, value_input, logic_select, spacing=5),
        )
        condition_btns = row(add_condition_btn, clear_conditions_btn, spacing=5)
        action_btns    = row(create_cluster_btn, delete_cluster_btn, spacing=5)

        _alias = getattr(self.component, "alias", None)
        _geom_label_div = Div(
            text="<b style='color:#888888'>Viewing:</b> "
                 + ((str(_alias) + "s") if _alias else "N/A")
        )

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
            _geom_label_div,
        ]

    # ------------------------------------------------------------------
    # Public: get_de_widgets (differential expression between groups)
    # ------------------------------------------------------------------

    def get_de_widgets(self, doc) -> list:
        self._de_groups_checkbox = CheckboxGroup(labels=[], active=[])
        self._de_groups_hint = Div(width=270)
        self._refresh_de_groups()

        def _on_de_groups_change(attr, old, new):
            names = self._de_groups_checkbox.labels
            self._de_selected = {names[i] for i in new}
        self._de_groups_checkbox.on_change("active", _on_de_groups_change)

        de_mode    = RadioButtonGroup(labels=["vs Rest", "Entre clusters"], active=0)
        run_de_btn = Button(label="▶ Run DE", button_type="primary", width=120)
        # Mensajes de estado/aviso del DE (Running DE..., bloqueos, sin resultados).
        self._de_status_div = Div(text="", width=270)

        def on_run_de():
            # In the order the groups were created: it decides overlaps in "Entre clusters"
            selected = [c for c in self._clusters if c in self._de_selected]
            if not selected:
                self._de_status_div.text = (
                    "<p style='color:red;font-size:11px'>"
                    "Select at least one group to compare.</p>"
                )
                return
            # "Entre clusters" (de_mode.active == 1) compara clusters entre sí:
            # con menos de 2 seleccionados no hay comparación posible y el hilo
            # de DE se quedaba colgado en "Running DE...". Bloquear aquí.
            if de_mode.active == 1 and len(selected) < 2:
                self._de_status_div.text = (
                    "<p style='color:red;font-size:11px'>"
                    "\"Entre clusters\" needs at least 2 selected groups.</p>"
                )
                return
            self._de_status_div.text = ""
            self._run_de_thread(selected, de_mode.active, doc)

        run_de_btn.on_click(on_run_de)

        return [
            Div(text="<b style='color:#888888'>Groups to compare</b>"),
            self._de_groups_hint,
            self._de_groups_checkbox,
            Div(text="<b style='color:#888888'>Compare against</b>"),
            de_mode,
            run_de_btn,
            self._de_status_div,
        ]

    def _refresh_de_groups(self):
        """Keep the DE group list in step with the groups that exist."""
        if self._de_groups_checkbox is None:
            return
        names = list(self._clusters)
        self._de_selected &= set(names)
        self._de_groups_checkbox.labels = names
        self._de_groups_checkbox.active = [i for i, c in enumerate(names) if c in self._de_selected]
        self._de_groups_hint.text = (
            "" if names else
            "<p style='color:#888;font-size:11px'>Create groups in the Groups tab first.</p>"
        )

    # ------------------------------------------------------------------
    # Cluster widgets (checkbox list)
    # ------------------------------------------------------------------

    def _build_cluster_widgets(self, doc) -> list:
        if not self._clusters:
            return [Div(
                text="<p style='color:#888;font-size:11px'>No clusters created yet.</p>"
            )]

        labels = []
        for cname, conditions in self._clusters.items():
            n_spots  = self._count_spots_bool(cname)
            # Cada condición (salvo la primera) va precedida de su lógica AND/OR.
            cond_str = ""
            for i, (col, op, val, logic) in enumerate(conditions):
                if i > 0:
                    cond_str += f" {logic} "
                cond_str += f"{col}{op}{val}"
            unit = f"{self.component.alias.lower()}s"  # spots / bins / cells
            labels.append(f"{cname}  ({n_spots} {unit})  [{cond_str}]")

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

    def _clear_de_if_cluster_deleted(self):
        """Quitar la tabla y los gráficos DE del plot si se borra alguno de los
        grupos con los que se calcularon (ocultarlos en el mapa no les afecta)."""
        de_clusters = getattr(self, '_de_clusters', None)
        if not de_clusters:
            return
        if all(c in self._clusters for c in de_clusters):
            return
        self._de_clusters = None
        self._de_results = {}
        self._last_plots_html = ""
        de_div = getattr(self.component, '_figure_de_div', None)
        if de_div is not None:
            de_div.text = ""
        if self._de_status_div is not None:
            self._de_status_div.text = ""

    def _rebuild_cluster_widgets(self, doc):
        if self._cluster_layout is None:
            return
        self._cluster_layout.children = self._build_cluster_widgets(doc)

    def _count_spots_bool(self, cname: str) -> int:
        col_key = f"grp_{cname}"
        obs = self.sample.adata.obs
        if col_key not in obs.columns:
            return 0
        return int(obs[col_key].astype(bool).sum())

    # ------------------------------------------------------------------
    # Visibility mask
    # ------------------------------------------------------------------

    def _apply_visibility_mask(self):
        if not self._clusters:
            return

        adata = self.sample.adata
        gdf   = self.sample.gdf
        mgr   = self.component.manager

        if hasattr(self.component, 'loader'):
            # Bins (lazy): el source solo tiene los bins cargados en el viewport,
            # así que la máscara debe alinearse con source.data['obs_name'] y no
            # con todo el gdf (si no, el BooleanFilter no casa y no se pinta nada).
            # Se usa la etiqueta de cluster que ya trae cada bin en el source: en
            # zoom alejado los bins agregados se llaman "summary:x:y:factor" y no
            # existen en adata.obs.
            n_source = len(mgr.source.data.get('x', []))
            labels = mgr.source.data.get(GROUPS_COL)
            if labels is not None and len(labels) == n_source:
                visible = pd.Series(labels, dtype=object).isin(self._active_clusters)
            else:
                obs_names = pd.Series(mgr.source.data.get('obs_name', []), dtype=object)
                visible = pd.Series(False, index=obs_names.index)
                for cname in self._active_clusters:
                    col_key = f"grp_{cname}"
                    if col_key in adata.obs.columns:
                        in_cluster = obs_names.map(adata.obs[col_key].astype(bool))
                        visible |= in_cluster.fillna(False).astype(bool)
            visible = visible.to_numpy(dtype=bool)
            if not self.component.show_not_sequenced:
                is_nan = mgr.get_nan_vector()
                if is_nan.shape[0] == visible.shape[0]:
                    visible &= ~is_nan
            mgr.view_filter = visible
        else:
            visible = pd.Series([False] * len(gdf), index=gdf.index)
            for cname in self._active_clusters:
                col_key = f"grp_{cname}"
                if col_key in adata.obs.columns:
                    cluster_series = self.sample.get_value(colname=col_key)
                    if not isinstance(cluster_series, pd.Series):
                        cluster_series = pd.Series(cluster_series, index=gdf.index)
                    cluster_series = cluster_series.reindex(gdf.index).fillna(False).astype(bool)
                    visible |= cluster_series

        col = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)
        for cname in self._clusters:
            if cname in self._active_clusters:
                col_key = f"grp_{cname}"
                if col_key in adata.obs.columns:
                    col[adata.obs[col_key].astype(bool)] = cname
        adata.obs[GROUPS_COL] = col

        self.component._view = visible.tolist()
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
                self._color_map[cname] = PALETTE[len(self._color_map) % len(PALETTE)]

        self._rebuild_custom_clusters_col()
        if hasattr(self.sample, 'clear_multiscale_cache'):
            self.sample.clear_multiscale_cache()
        try:
            self.component.clustering_manager.apply_colors_to_column(GROUPS_COL)
        except ValueError:
            color_param = self.component.manager.get_param('color')
            if color_param is not None:
                color_param.colname = GROUPS_COL
                self.component.manager.update_param('color')

        color_param = self.component.manager.get_param('color')
        if color_param is not None and hasattr(color_param, 'category_colors'):
            color_param.category_colors = dict(self._color_map)
            if hasattr(self.component, 'loader'):
                self.component.manager.update_param('alpha', avoid_scaling=True)
                self.component.loader.force_reload()
                return
            obs_series = self.sample.adata.obs[GROUPS_COL]
            obs_names  = pd.Series(self.component.manager.source.data['obs_name'])
            aligned    = obs_names.map(obs_series)
            self.component.manager.cache[GROUPS_COL] = aligned.values
            new_colors = color_param.get_values(
                pd.Series(self.component.manager.cache[GROUPS_COL]).copy()
            )
            data = dict(self.component.manager.source.data)
            data['color']           = new_colors
            data['custom_clusters'] = ['—' if pd.isna(c) else c for c in aligned]
            self.component.manager.source.data = data
            self.component.manager.fetch_data()

    def _rebuild_custom_clusters_col(self):
        adata = self.sample.adata
        col   = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)
        for cname in self._clusters:
            col_key = f"grp_{cname}"
            if col_key in adata.obs.columns:
                col[adata.obs[col_key].astype(bool)] = cname
        adata.obs[GROUPS_COL] = col.astype('category')

    # ------------------------------------------------------------------
    # Build boolean mask from conditions
    # ------------------------------------------------------------------

    def _build_mask(self, conditions, adata) -> pd.Series:
        obs     = adata.obs
        var_idx = {g: i for i, g in enumerate(adata.var_names)}
        mask    = None

        for col, op, val, logic in conditions:
            val_f = float(val)

            if col in var_idx:
                expr = adata.X[:, var_idx[col]]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                series = pd.Series(expr, index=obs.index)
            elif col in obs.columns:
                series = obs[col].astype(float)
            else:
                continue

            ops       = {'>': series > val_f, '<': series < val_f,
                         '>=': series >= val_f, '<=': series <= val_f}
            cond_mask = ops.get(op, series == val_f)

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
