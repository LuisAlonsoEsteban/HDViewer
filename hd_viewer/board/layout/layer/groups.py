
from html import escape

import numpy as np
import pandas as pd

from bokeh.models import (
    Button, Div, Select, TextInput,
    AutocompleteInput, CheckboxGroup, ColorPicker,
)
from bokeh.layouts import column, row

from .... import settings
from ...utils.loading import LoadingLogger
from ....startup.utils.progress import report, reporting

from ....core.groups import GROUPS_COL, OPERATORS, GroupsState  # noqa: F401  (re-exported)

class GroupsManager(GroupsState):
    """
    Manages creation of custom groups based on gene/metadata conditions.
    Stores results in adata.obs[GROUPS_COL], same structure as leiden.

    Each unit (bins/cells) gets the name of the first matching group, or None if unassigned.
    """

    def __init__(self, component):
        self.component = component
        super().__init__()
        self._conditions: list = []   # conditions being built [(col, op, val, logic)]

        self._group_checkbox = None
        self._group_layout = None

        self._conditions_div = None
        self._status_div = None
        self.logger = None
        self._operation_running = False

    @property
    def _base_sample(self):
        return self.component.manager.sample






    def add_mask_group(self, name, mask, doc, origin_roi_id=None, prepared=None):
        """Create a fixed membership group, independently of gene conditions."""
        existing = self.group_for_roi(origin_roi_id) if origin_roi_id is not None else None
        name = self.add_fixed_group(name, mask, origin_roi_id=origin_roi_id, prepared=prepared)
        if existing is None:
            self._rebuild_group_widgets(doc)
            self._notify_group_change()
        return name

    def delete_group(self, name, doc, prepared=None):
        """Shared deletion path for condition groups and ROI snapshot groups."""
        if name not in self._groups:
            return
        members = self.members(name)
        if members:
            self.component.set_selected_obs_names(self.component.get_selected_obs_names() - members)
        if self.forget_group(name):
            self._apply_groups_to_adata(prepared)
        else:
            if GROUPS_COL in self.sample.adata.obs:
                self.component.show_colname(self.sample.adata.var_names[0])
                self.sample.adata.obs.drop(columns=[GROUPS_COL], inplace=True)
            manager = self.component.manager
            visible = np.ones(len(manager.loaded_obs_names), dtype=bool)
            if not self.component.show_not_sequenced:
                visible &= ~manager.get_nan_vector()
            if self.component.lazy_loading:
                manager.view_filter = visible
            self.component._view = visible
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

        self.logger = LoadingLogger(doc, width=270, font_size=11)
        self._status_div = self.logger.div

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
            if self._operation_running:
                return
            col = self.gene_metadata_input.value.strip()
            op = operator_select.value
            val = value_input.value.strip()
            logic = logic_select.value if self._conditions else "AND"

            if not col or not val:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Gene and value are required.</p>"
                )
                return


            try:
                hint = self.check_condition(col, op, val)
            except ValueError as exc:
                self._status_div.text = f"<p style='color:red;font-size:11px'>{escape(str(exc))}</p>"
                return
            self._status_div.text = f"<p style='color:#4e9af1;font-size:11px'>{escape(hint)}</p>" if hint else ""

            self._conditions.append((col, op, val, logic))
            doc.add_next_tick_callback(_render_conditions)

        def on_clear_conditions():
            if self._operation_running:
                return
            self._conditions.clear()
            doc.add_next_tick_callback(_render_conditions)

        def on_create_group():
            if self._operation_running:
                return
            if not self._conditions:
                self._status_div.text = (
                    "<p style='color:red;font-size:11px'>Add at least one condition first.</p>"
                )
                return

            cname = group_name_input.value.strip() or "Group A"
            if cname in self._groups:
                self._status_div.text = (
                    f"<p style='color:red;font-size:11px'>Group '{escape(cname)}' already exists.</p>"
                )
                return

            if self.logger.running:
                return
            conditions = list(self._conditions)

            def done(result):
                self.commit_condition_group(*result)
                self._conditions.clear()
                _render_conditions()
                self._rebuild_group_widgets(doc)

            self.run_group_action(doc, f"Creating group '{cname}'",
                                  lambda: self.prepare_condition_group(cname, conditions), done,
                                  f"Group '{escape(cname)}' created.")

        def on_delete_group():
            cname = group_name_input.value.strip()
            if self.logger.running:
                return
            if cname not in self._groups:
                self.logger.error(f"Group '{cname}' does not exist.")
                return
            remaining = {name: conditions for name, conditions in self._groups.items() if name != cname}
            self.run_group_action(doc, f"Removing group '{cname}'",
                                  lambda: self._prepare_groups(remaining),
                                  lambda prepared: self.delete_group(cname, doc, prepared),
                                  f"Group '{escape(cname)}' deleted.")

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
        color_pickers = []
        for name in self._groups:
            picker = ColorPicker(title=name, color=self._color_map[name], width=130)

            def on_color_change(attr, old, new, group_name=name):
                if group_name not in self._groups:
                    return
                self._color_map[group_name] = new
                manager = self.component.manager
                param = manager.get_param('color')
                if getattr(param, 'colname', None) == GROUPS_COL:
                    param.category_colors = dict(self._color_map)
                    data = dict(manager.source.data)
                    if self.component.lazy_loading and GROUPS_COL in data:
                        labels = pd.Series(data[GROUPS_COL])
                    else:
                        labels = pd.Series(manager.loaded_obs_names, dtype=object).map(
                            self.sample.adata.obs[GROUPS_COL]
                        )
                    data['color'] = param.get_values(labels)
                    manager.source.data = data

            picker.on_change('color', on_color_change)
            color_pickers.append(picker)
        return [self._group_checkbox, *color_pickers]

    def _rebuild_group_widgets(self, doc):
        if self._group_layout is None:
            return
        self._group_layout.children = self._build_group_widgets(doc)


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
                visible = pd.Series(self.visible_mask(mgr.loaded_obs_names))
            visible = visible.to_numpy(dtype=bool)
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

        # Unsequenced geometry is independent of group membership.
        visible = np.asarray(visible, dtype=bool)
        is_nan = mgr.get_nan_vector()
        visible[is_nan] = self.component.show_not_sequenced
        if self.component.lazy_loading:
            mgr.view_filter = visible

        self.refresh_active_labels()

        self.component._view = visible
        self.component.change_glyph_view()

    # ------------------------------------------------------------------
    # Apply groups to adata.obs[GROUPS_COL]
    # ------------------------------------------------------------------

    def run_group_action(self, doc, title, work, done, success, logger=None):
        """Compute membership in a worker; update Bokeh models on the document thread."""
        logger = logger or self.logger
        if logger is None:
            self.logger = logger = LoadingLogger(doc, width=270, font_size=11)
        if logger.running or self._operation_running:
            return

        self._operation_running = True

        def failed(exc):
            self._operation_running = False
            logger.error(exc)

        def finish(result):
            try:
                with reporting(logger.log):
                    report("Updating group colours and map...")
                    done(result)
                logger.success(success)
            except Exception as exc:
                logger.error(exc)
            finally:
                self._operation_running = False

        logger.run(title, work, finish, failed)


    def _apply_groups_to_adata(self, prepared=None):
        """Store the groups, then colour the map by them."""
        self.store_groups(prepared)
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
            obs_names  = pd.Series(self.component.manager.loaded_obs_names)
            aligned    = obs_names.map(obs_series)
            self.component.manager.cache[GROUPS_COL] = aligned.values
            new_colors = color_param.get_values(
                pd.Series(self.component.manager.cache[GROUPS_COL]).copy()
            )
            data = dict(self.component.manager.source.data)
            data['color']           = new_colors
            data['custom_groups'] = aligned.astype(object).where(aligned.notna(), '—').tolist()
            self.component.manager.source.data = data
            self.component.manager.fetch_data()


    # ------------------------------------------------------------------
    # Build boolean mask from conditions
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
