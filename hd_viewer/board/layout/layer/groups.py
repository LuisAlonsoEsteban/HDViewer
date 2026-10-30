
import numpy as np
import pandas as pd

from bokeh.models import (
    Button, Div, Select, TextInput,
    AutocompleteInput, CheckboxGroup,
)
from bokeh.layouts import column, row

from .... import settings
from ...utils.palette import GROUP_PALETTE, palette_color

OPERATORS  = [">", "<", ">=", "<=", "=="]
GROUPS_COL = settings.GROUPS_COL

class GroupsManager:
    """
    Manages creation of custom groups based on gene/metadata conditions.
    Stores results in adata.obs[GROUPS_COL], same structure as leiden.

    Each unit (bins/cells) gets the name of the first matching group, or None if unassigned.
    """

    def __init__(self, component):
        self.component = component

        self._groups: dict = {}   # group_name → [conditions]
        self._active_groups: set = set()

        self._fixed_masks = {}
        self._group_origins = {}
        self._change_listeners = []
        self._conditions: list = []   # conditions being built [(col, op, val, logic)]
        
        self._group_checkbox = None
        self._group_layout = None

        self._color_map: dict = {}

        self._conditions_div = None
        self._status_div = None

    @property
    def sample(self):
        return self.component.manager.sample

    @property
    def has_groups(self) -> bool:
        return (
            GROUPS_COL in self.sample.adata.obs.columns
            and len(self._groups) > 0
        )

    def add_change_listener(self, callback):
        if callback not in self._change_listeners:
            self._change_listeners.append(callback)

    def remove_change_listener(self, callback):
        if callback in self._change_listeners:
            self._change_listeners.remove(callback)

    def _notify_group_change(self):
        for callback in list(self._change_listeners):
            callback()

    def group_for_roi(self, roi_id):
        return next((name for name, origin in self._group_origins.items() if origin == roi_id), None)

    def add_mask_group(self, name, mask, doc, origin_roi_id=None):
        """Create a fixed membership group, independently of gene conditions."""
        mask = np.asarray(mask, dtype=bool)
        
        if mask.shape != (self.sample.adata.n_obs,) or not mask.any():
            raise ValueError("A group needs a nonempty mask aligned with the sample observations.")
        if origin_roi_id is not None and self.group_for_roi(origin_roi_id):
            return self.group_for_roi(origin_roi_id)
        
        base, suffix = name, 2
        while name in self._groups:
            name = f"{base} ({suffix})"
            suffix += 1

        self._fixed_masks[name] = pd.Series(mask.copy(), index=self.sample.adata.obs_names.copy())
        self._groups[name] = []
        self._group_origins[name] = origin_roi_id
        self._active_groups.add(name)
        
        try:
            self._apply_groups_to_adata()
        except Exception:
            self._groups.pop(name)
            self._fixed_masks.pop(name)
            self._group_origins.pop(name)
            self._active_groups.discard(name)
            self.sample.adata.obs.drop(columns=[f"{settings.GROUP_MASK_PREFIX}{name}"], errors="ignore", inplace=True)
            raise
        
        self._rebuild_group_widgets(doc)
        self._notify_group_change()
        return name

    def delete_group(self, name, doc):
        """Shared deletion path for condition groups and ROI snapshot groups."""
        if name not in self._groups:
            return
        del self._groups[name]
        self._fixed_masks.pop(name, None)
        self._group_origins.pop(name, None)
        self._active_groups.discard(name)
        self.sample.adata.obs.drop(columns=[f"{settings.GROUP_MASK_PREFIX}{name}"], errors="ignore", inplace=True)
        if hasattr(self.sample, 'clear_multiscale_cache'):
            self.sample.clear_multiscale_cache()
        if self._groups:
            self._apply_groups_to_adata()
        else:
            if GROUPS_COL in self.sample.adata.obs:
                self.component.show_colname(self.sample.adata.var_names[0])
                self.sample.adata.obs.drop(columns=[GROUPS_COL], inplace=True)
            manager = self.component.manager
            visible = np.ones(len(manager.source.data.get('obs_name', [])), dtype=bool)
            if not self.component.show_not_sequenced:
                visible &= ~manager.get_nan_vector()
            if self.component.lazy_loading:
                manager.view_filter = visible
            self.component._view = visible.tolist()
            self.component.change_glyph_view()
            if hasattr(self.component, 'loader'):
                self.component.loader.force_reload()
            else:
                manager.fetch_data()
        self._rebuild_group_widgets(doc)
        self._notify_group_change()

    def refresh_completions(self):
        if getattr(self, "gene_metadata_input", None):
            self.gene_metadata_input.completions = self.sample.get_colnames()

    def get_widgets(self, doc) -> list:
        colnames = self.sample.get_colnames()
        self.gene_metadata_input = AutocompleteInput(
            title="Gene / Metadata",
            value=colnames[0] if colnames else "",
            completions=colnames,
            min_characters=settings.METADATA_INPUT_MIN_CHARACTERS,
            case_sensitive=settings.AUTOCOMPLETE_CASE_SENSITIVE,
            search_strategy=settings.AUTOCOMPLETE_SEARCH_STRATEGY,
            restrict=False,
            width=235,
        )
        operator_select = Select(title="Operator", options=OPERATORS, value=">", width=70)
        value_input = TextInput(title="Value", value="1", width=80)
        logic_select = Select(title="Logic", options=["AND", "OR"], value="AND", width=70)

        self._conditions_div = Div(
            text="<p style='color:#888;font-size:11px'>No conditions added yet.</p>",
            width=270,
        )
        group_name_input = TextInput(value="Group A", title="Group name", width=200)

        add_condition_btn = Button(label="+ Add condition",  button_type="success", width=150)
        clear_conditions_btn = Button(label="Clear", button_type="warning",  width=80)
        create_group_btn = Button(label="✓ Create group", button_type="success", width=130)
        delete_group_btn = Button(label="🗑 Delete group", button_type="danger",  width=130)

        self._status_div = Div(text="", width=270)

        if self._group_layout is None:
            self._group_layout = column(
                name="groups_controls",
                children=self._build_group_widgets(doc),
            )

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

        def on_add_condition():
            col = self.gene_metadata_input.value.strip()
            op = operator_select.value
            val = value_input.value.strip()
            logic = logic_select.value if self._conditions else "AND"

            if not col or not val:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Gene and value are required.</p>"
                )
                return


            if col in self.sample.adata.var_names:
                try:
                    float(val)
                except ValueError:
                    self._status_div.text = (
                        "<p style='color:red;font-size:11px'>Gene expression conditions require a numeric value.</p>"
                    )
                    return

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
            elif col in self.sample.adata.obs.columns:
                series = self.sample.adata.obs[col]
                if pd.api.types.is_numeric_dtype(series):
                    try:
                        float(val)
                    except ValueError:
                        self._status_div.text = (
                            "<p style='color:red;font-size:11px'>Numeric metadata requires a numeric value.</p>"
                        )
                        return
                elif op in (">", "<", ">=", "<="):
                    self._status_div.text = (
                        "<p style='color:red;font-size:11px'>Only == is supported for non-numeric metadata values.</p>"
                    )
                    return
                self._status_div.text = ""
            else:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Unknown gene or metadata column.</p>"
                )
                return

            self._conditions.append((col, op, val, logic))
            doc.add_next_tick_callback(_render_conditions)

        def on_clear_conditions():
            self._conditions.clear()
            doc.add_next_tick_callback(_render_conditions)

        def on_create_group():
            if not self._conditions:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Add at least one condition first.</p>"
                )
                return

            cname = group_name_input.value.strip() or "Group A"
            if cname in self._groups:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Group '{cname}' already exists.</p>"
                )
                return

            self._groups[cname] = list(self._conditions)
            self._active_groups.add(cname)
            self._conditions.clear()

            try:
                self._apply_groups_to_adata()
                self._status_div.text = (
                    f"<p style='color:green;font-size:11px'>Group '{cname}' created.</p>"
                )
            except Exception as e:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Error: {str(e)[:80]}</p>"
                )

            def _update():
                _render_conditions()
                self._rebuild_group_widgets(doc)
                self._notify_group_change()

            doc.add_next_tick_callback(_update)

        def on_delete_group():
            cname = group_name_input.value.strip()
            if cname not in self._groups:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Group '{cname}' does not exist.</p>"
                )
                return

            try:
                self.delete_group(cname, doc)
                self._status_div.text = f"<p style='color:#f0a500'>Group '{cname}' deleted.</p>"
            except Exception as exc:
                self._status_div.text = f"<p style='color:red'>Error: {str(exc)[:80]}</p>"

        add_condition_btn.on_click(on_add_condition)
        clear_conditions_btn.on_click(on_clear_conditions)
        create_group_btn.on_click(on_create_group)
        delete_group_btn.on_click(on_delete_group)

        condition_row  = column(
            self.gene_metadata_input,
            row(operator_select, value_input, logic_select, spacing=5),
        )
        condition_btns = row(add_condition_btn, clear_conditions_btn, spacing=5)
        action_btns    = row(create_group_btn, delete_group_btn, spacing=5)

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
            group_name_input,
            action_btns,
            self._status_div,
            Div(text="<b style='color:#888888'>Created groups</b>"),
            self._group_layout,
            _geom_label_div,
        ]

    # ------------------------------------------------------------------
    # Group widgets (checkbox list)
    # ------------------------------------------------------------------

    def _build_group_widgets(self, doc) -> list:
        if not self._groups:
            return [Div(
                text="<p style='color:#888;font-size:11px'>No groups created yet.</p>"
            )]

        labels = []
        for cname, conditions in self._groups.items():
            n_spots  = self._count_spots_bool(cname)
            cond_str = ""
            for i, (col, op, val, logic) in enumerate(conditions):
                if i > 0:
                    cond_str += f" {logic} "
                cond_str += f"{col}{op}{val}"
            if cname in self._fixed_masks:
                cond_str = "ROI membership"
            unit = f"{self.component.alias.lower()}s"
            labels.append(f"{cname}  ({n_spots} {unit})  [{cond_str}]")

        self._group_checkbox = CheckboxGroup(
            labels=labels,
            active=[i for i, cname in enumerate(self._groups)
                    if cname in self._active_groups],
        )
        self._apply_visibility_mask()

        def _on_checkbox_change(attr, old, new):
            group_names = list(self._groups.keys())
            self._active_groups = {group_names[i] for i in new}
            doc.add_next_tick_callback(self._apply_visibility_mask)

        self._group_checkbox.on_change("active", _on_checkbox_change)
        return [self._group_checkbox]

    def _rebuild_group_widgets(self, doc):
        if self._group_layout is None:
            return
        self._group_layout.children = self._build_group_widgets(doc)

    def _count_spots_bool(self, cname: str) -> int:
        col_key = f"{settings.GROUP_MASK_PREFIX}{cname}"
        obs = self.sample.adata.obs
        if col_key not in obs.columns:
            return 0
        return int(obs[col_key].astype(bool).sum())

    # ------------------------------------------------------------------
    # Visibility mask
    # ------------------------------------------------------------------

    def _apply_visibility_mask(self):
        if not self._groups:
            return

        adata = self.sample.adata
        gdf   = self.sample.gdf
        mgr   = self.component.manager

        if hasattr(self.component, 'loader'):
            n_source = len(mgr.source.data.get('x', []))
            labels = mgr.source.data.get(GROUPS_COL)
            if labels is not None and len(labels) == n_source:
                visible = pd.Series(labels, dtype=object).isin(self._active_groups)
            else:
                obs_names = pd.Series(mgr.source.data.get('obs_name', []), dtype=object)
                visible = pd.Series(False, index=obs_names.index)
                for gname in self._active_groups:
                    col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
                    if col_key in adata.obs.columns:
                        in_group = obs_names.map(adata.obs[col_key].astype(bool))
                        visible |= in_group.fillna(False).astype(bool)
            visible = visible.to_numpy(dtype=bool)
            if not self.component.show_not_sequenced:
                is_nan = mgr.get_nan_vector()
                if is_nan.shape[0] == visible.shape[0]:
                    visible &= ~is_nan
            mgr.view_filter = visible
        else:
            visible = pd.Series([False] * len(gdf), index=gdf.index)
            for gname in self._active_groups:
                col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
                if col_key in adata.obs.columns:
                    group_series = self.sample.get_value(colname=col_key)
                    if not isinstance(group_series, pd.Series):
                        group_series = pd.Series(group_series, index=gdf.index)
                    group_series = group_series.reindex(gdf.index).fillna(False).astype(bool)
                    visible |= group_series

        col = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)
        for gname in self._groups:
            if gname in self._active_groups:
                col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
                if col_key in adata.obs.columns:
                    col[adata.obs[col_key].astype(bool)] = gname
        adata.obs[GROUPS_COL] = col

        self.component._view = visible.tolist()
        self.component.change_glyph_view()

    # ------------------------------------------------------------------
    # Apply groups to adata.obs[GROUPS_COL]
    # ------------------------------------------------------------------

    def _apply_groups_to_adata(self):
        adata = self.sample.adata

        for gname, conditions in self._groups.items():
            mask = (self._fixed_masks[gname].reindex(adata.obs_names, fill_value=False)
                    if gname in self._fixed_masks else self._build_mask(conditions, adata))
            adata.obs[f"{settings.GROUP_MASK_PREFIX}{gname}"] = mask.values
            if gname not in self._color_map:
                self._color_map[gname] = palette_color(GROUP_PALETTE, len(self._color_map))

        self._rebuild_custom_groups_col()
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
            data['custom_groups'] = ['—' if pd.isna(c) else c for c in aligned]
            self.component.manager.source.data = data
            self.component.manager.fetch_data()

    def _rebuild_custom_groups_col(self):
        adata = self.sample.adata
        col   = pd.Series([None] * len(adata.obs), index=adata.obs.index, dtype=object)
        for gname in self._groups:
            col_key = f"{settings.GROUP_MASK_PREFIX}{gname}"
            if col_key in adata.obs.columns:
                col[adata.obs[col_key].astype(bool)] = gname
        adata.obs[GROUPS_COL] = col.astype('category')

    # ------------------------------------------------------------------
    # Build boolean mask from conditions
    # ------------------------------------------------------------------

    def _build_mask(self, conditions, adata) -> pd.Series:
        obs     = adata.obs
        var_idx = {g: i for i, g in enumerate(adata.var_names)}
        mask    = None

        for col, op, val, logic in conditions:
            if col in var_idx:
                expr = adata.X[:, var_idx[col]]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                series = pd.Series(expr, index=obs.index)
                try:
                    val_f = float(val)
                except ValueError as exc:
                    raise ValueError(f"Gene condition '{col}' requires a numeric value.") from exc
                ops = {'>': series > val_f, '<': series < val_f,
                       '>=': series >= val_f, '<=': series <= val_f}
                cond_mask = ops.get(op, series == val_f)
            elif col in obs.columns:
                series = obs[col]
                if pd.api.types.is_numeric_dtype(series):
                    try:
                        val_f = float(val)
                    except ValueError as exc:
                        raise ValueError(f"Column '{col}' requires a numeric value for operator '{op}'.") from exc
                    ops = {'>': series > val_f, '<': series < val_f,
                           '>=': series >= val_f, '<=': series <= val_f}
                    cond_mask = ops.get(op, series == val_f)
                else:
                    if op not in ("==",):
                        raise ValueError(f"Operator '{op}' is only supported for numeric columns; use == for '{col}'.")
                    cond_mask = series.astype(str) == str(val)
            else:
                continue

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
