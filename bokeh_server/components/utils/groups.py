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
            width=160,
        )
        operator_select = Select(title="Operator", options=OPERATORS, value=">", width=70)
        value_input     = TextInput(title="Value", value="1", width=80)
        logic_select    = Select(title="Logic", options=["AND", "OR"], value="AND", width=70)

        self._conditions_div = Div(
            text="<p style='color:#888;font-size:11px'>No conditions added yet.</p>",
            width=380,
        )
        cluster_name_input = TextInput(value="Cluster A", title="Cluster name", width=200)

        add_condition_btn    = Button(label="+ Add condition",  button_type="success", width=150)
        clear_conditions_btn = Button(label="Clear",            button_type="warning",  width=80)
        create_cluster_btn   = Button(label="✓ Create cluster", button_type="success", width=160)
        delete_cluster_btn   = Button(label="🗑 Delete cluster", button_type="danger",  width=160)
        de_mode              = RadioButtonGroup(labels=["vs Rest", "Entre clusters"], active=0)
        run_de_btn           = Button(label="▶ Run DE", button_type="primary", width=120)

        self._status_div = Div(text="", width=380)

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
            self._conditions.clear()

            try:
                self._apply_clusters_to_adata()
                self._status_div.text = (
                    f"<p style='color:green;font-size:11px'>Cluster '{cname}' created.</p>"
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
                    f"<p style='color:red;font-size:11px'>Cluster '{cname}' does not exist.</p>"
                )
                return

            del self._clusters[cname]
            self._active_clusters.discard(cname)

            col_key = f"grp_{cname}"
            if col_key in self.sample.adata.obs.columns:
                self.sample.adata.obs.drop(columns=[col_key], inplace=True)

            try:
                if self._clusters:
                    self._apply_clusters_to_adata()
                else:
                    if GROUPS_COL in self.sample.adata.obs.columns:
                        self.sample.adata.obs.drop(columns=[GROUPS_COL], inplace=True)
                self._status_div.text = (
                    f"<p style='color:#f0a500;font-size:11px'>Cluster '{cname}' deleted.</p>"
                )
            except Exception as e:
                import traceback
                traceback.print_exc()
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            doc.add_next_tick_callback(lambda: self._rebuild_cluster_widgets(doc))

        def on_run_de():
            active = list(self._active_clusters)
            if active:
                self._run_de_thread(active, de_mode.active, doc)

        add_condition_btn.on_click(on_add_condition)
        clear_conditions_btn.on_click(on_clear_conditions)
        create_cluster_btn.on_click(on_create_cluster)
        delete_cluster_btn.on_click(on_delete_cluster)
        run_de_btn.on_click(on_run_de)

        # --- Layout ---
        condition_row  = row(gene_input, operator_select, value_input, logic_select, spacing=5)
        condition_btns = row(add_condition_btn, clear_conditions_btn, spacing=5)
        action_btns    = row(create_cluster_btn, delete_cluster_btn, spacing=10)

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
            Div(text="<hr style='border-color:#444'>"),
            de_mode,
            run_de_btn,
            _geom_label_div,
        ]

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
            cond_str = " & ".join(f"{col}{op}{val}" for col, op, val, logic in conditions)
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
                print(f"[GroupsManager] Warning: '{col}' not found, skipping.")
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
