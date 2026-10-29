
from html import escape

import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import (
    AutocompleteInput,
    Button,
    CheckboxGroup,
    ColorPicker,
    Div,
    NumericInput,
    RadioButtonGroup,
    Select,
)

from .... import settings
from ....help import help_header
from ...source.parameter import CategoricalColor, Scaler, StaticParameter


CORE_ATTRS = ("color", "alpha")
MODES = ("static", "data")

MUTED = "color:#888;font-size:11px"

def is_color(attr):
    return "color" in attr


def attr_spec(attr):
    """Input limits and defaults per kind of attribute."""
    if is_color(attr):
        return dict(settings.PARAMETER_COLOR_SPEC)
    if "alpha" in attr:
        return dict(settings.PARAMETER_ALPHA_SPEC)
    return dict(settings.PARAMETER_WIDTH_SPEC)  # widths


class ParamForm:

    def __init__(self, component):
        self.component = component
        self._editing = None   # attr being edited, or None when adding

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def get_layout(self, doc):
        self.list_column = column(sizing_mode="stretch_width", spacing=4)
        self.form_title = Div()

        self.attr_select = Select(title="Attribute", width=235)
        self.mode = RadioButtonGroup(labels=["Fixed value", "From data"], active=0)

        # Fixed value
        self.value_color = ColorPicker(title="Value", width=235)
        self.value_number = NumericInput(title="Value", mode="float", width=235)
        # From data
        self.colname_input = AutocompleteInput(
            title="Gene / Metadata", value="", min_characters=settings.METADATA_INPUT_MIN_CHARACTERS,
            case_sensitive=settings.AUTOCOMPLETE_CASE_SENSITIVE,
            search_strategy=settings.AUTOCOMPLETE_SEARCH_STRATEGY, restrict=True, width=235,
        )
        self.colname_hint = Div(width=235)
        self.min_input = NumericInput(title="Output min", mode="float", width=110)
        self.max_input = NumericInput(title="Output max", mode="float", width=110)
        self.range_row = row(self.min_input, self.max_input, spacing=10)
        self.hover_checkbox = CheckboxGroup(labels=["Show the column in the hover"], active=[])
        # Both modes
        self.ns_color = ColorPicker(title="Value for not-sequenced", width=235)
        self.ns_number = NumericInput(title="Value for not-sequenced", mode="float", width=235)

        self.submit_button = Button(label="＋ Add parameter", button_type="success", width=150)
        self.cancel_button = Button(label="Cancel", width=80, visible=False)
        self.status = Div(width=270)

        self.attr_select.on_change("value", lambda attr, old, new: self._sync(reset_values=True))
        self.mode.on_change("active", lambda attr, old, new: self._sync(reset_values=False))
        self.submit_button.on_click(self._on_submit)
        self.cancel_button.on_click(self._start_adding)

        self.form = column(
            self.form_title, self.attr_select, self.mode,
            self.value_color, self.value_number,
            self.colname_input, self.colname_hint, self.range_row, self.hover_checkbox,
            self.ns_color, self.ns_number,
            row(self.submit_button, self.cancel_button, spacing=6),
            self.status,
            spacing=6,
        )
        self.layout = column(
            help_header("Parameters", "parameters"),
            Div(text=f"<p style='{MUTED};margin:0'>Colour and alpha are set in the "
                     "Color and General tabs.</p>"),
            self.list_column,
            Div(text="<hr style='border-color:#444'>"),
            self.form,
            spacing=6,
        )
        self._start_adding()
        return self.layout

    # ------------------------------------------------------------------
    # Current parameters
    # ------------------------------------------------------------------

    def extra_parameters(self):
        return [p for p in self.component.manager.parameters if p.attr not in CORE_ATTRS]

    @staticmethod
    def describe(param):
        def swatch(value):
            return (f"<span style='display:inline-block;width:10px;height:10px;"
                    f"background:{escape(str(value))};border:1px solid #888;"
                    f"vertical-align:middle;margin-right:3px'></span>")

        def fmt(value, attr):
            return swatch(value) + escape(str(value)) if is_color(attr) else f"{value:g}"

        if isinstance(param, StaticParameter):
            how = f"= {fmt(param.value, param.attr)}"
        elif isinstance(param, CategoricalColor):
            how = f"← {escape(param.colname)} <span style='color:#888'>(one colour per category)</span>"
        else:
            how = (f"← {escape(param.colname)} "
                   f"<span style='color:#888'>(scaled {param.minv:g}–{param.maxv:g})</span>")
        ns = fmt(param.not_sequenced_val, param.attr)
        return (f"<b>{escape(param.attr)}</b> {how}<br>"
                f"<span style='{MUTED}'>not sequenced: {ns}</span>")

    def _refresh_list(self):
        rows = []
        for param in self.extra_parameters():
            edit = Button(label="✎", width=32, height=26, margin=(0, 2, 0, 0))
            remove = Button(label="✕", width=32, height=26, button_type="danger", margin=0)
            edit.on_click(lambda attr=param.attr: self._start_editing(attr))
            remove.on_click(lambda attr=param.attr: self._remove(attr))
            rows.append(row(Div(text=self.describe(param), width=200,
                                styles=dict(color="#e5edf0", font_size="12px")),
                            edit, remove, spacing=4))
        if not rows:
            rows = [Div(text=f"<p style='{MUTED};margin:0'>No extra parameters yet.</p>")]
        self.list_column.children = rows

    def _available_attrs(self):
        used = {p.attr for p in self.component.manager.parameters}
        return [a for a in self.component.AVAILABLE_ATTRS if a not in CORE_ATTRS and a not in used]

    # ------------------------------------------------------------------
    # Form state
    # ------------------------------------------------------------------

    def _start_adding(self):
        self._editing = None
        options = self._available_attrs()
        self.attr_select.disabled = False
        self.attr_select.options = options
        self.form_title.text = "<b style='color:#888888'>Add a parameter</b>"
        self.submit_button.label = "＋ Add parameter"
        self.cancel_button.visible = False
        self.form.visible = bool(options)
        self._refresh_list()
        if options:
            self.mode.active = 0
            if self.attr_select.value in options:
                self._sync(reset_values=True)
            else:
                self.attr_select.value = options[0]   # triggers _sync

    def _start_editing(self, attr):
        param = self.component.manager.get_param(attr)
        if param is None:
            return
        self._editing = attr
        self.form.visible = True
        self.form_title.text = f"<b style='color:#888888'>Edit {escape(attr)}</b>"
        self.submit_button.label = "Save changes"
        self.cancel_button.visible = True
        self.status.text = ""
        self.attr_select.options = [attr]
        self.attr_select.value = attr
        self.attr_select.disabled = True
        self.mode.active = 0 if param.is_static else 1
        self._sync(reset_values=False)

        if param.is_static:
            self._set_value(param.value)
        else:
            self.colname_input.value = param.colname
            self.hover_checkbox.active = [0] if param.colname_hovered else []
            if isinstance(param, Scaler):
                self.min_input.value, self.max_input.value = param.minv, param.maxv
        self._set_not_sequenced(param.not_sequenced_val)

    def _sync(self, reset_values):
        """Show only the inputs that apply to the attribute and the mode."""
        attr = self.attr_select.value
        if not attr:
            return
        color, static = is_color(attr), self.mode.active == 0
        spec = attr_spec(attr)

        self.value_color.visible = static and color
        self.value_number.visible = static and not color
        self.colname_input.visible = not static
        self.colname_hint.visible = not static
        self.range_row.visible = not static and not color
        self.hover_checkbox.visible = not static
        self.ns_color.visible = color
        self.ns_number.visible = not color

        for widget in (self.value_number, self.ns_number, self.min_input, self.max_input):
            widget.low, widget.high = spec.get("low"), spec.get("high")
        self.colname_input.completions = self._column_choices(attr)
        self.colname_hint.text = (
            f"<p style='{MUTED};margin:0'>"
            + ("Categorical metadata only: each category gets a colour."
               if color else "A gene or numeric metadata, scaled to the output range.")
            + "</p>")
        if reset_values:
            self._set_value(spec["value"])
            self._set_not_sequenced(spec["not_sequenced"])
            self.colname_input.value = ""
            self.hover_checkbox.active = []
            self.min_input.value, self.max_input.value = spec.get("minv"), spec.get("maxv")
            self.status.text = ""

    def _column_choices(self, attr):
        adata = self.component.sample.adata
        numeric = {c: pd.api.types.is_numeric_dtype(adata.obs[c]) for c in adata.obs.columns}
        if is_color(attr):
            return [c for c, is_num in numeric.items() if not is_num]
        # Genes and numeric metadata (what Sample.is_continous accepts)
        return list(adata.var_names) + [c for c, is_num in numeric.items() if is_num]

    def refresh_completions(self):
        """New columns appeared (e.g. clustering): offer them."""
        if self.colname_input and self.attr_select.value:
            self.colname_input.completions = self._column_choices(self.attr_select.value)

    def _set_value(self, value):
        if is_color(self.attr_select.value):
            self.value_color.color = value
        else:
            self.value_number.value = value

    def _set_not_sequenced(self, value):
        if is_color(self.attr_select.value):
            self.ns_color.color = value
        else:
            self.ns_number.value = value

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def _error(self, message):
        self.status.text = f"<p style='color:#fb9292;font-size:11px;margin:0'>{escape(message)}</p>"

    def retrieve_form_values(self):
        """Arguments for Component.add_parameter; raises ValueError on invalid input."""
        attr = self.attr_select.value
        if not attr:
            raise ValueError("Choose an attribute.")
        color, static = is_color(attr), self.mode.active == 0
        values = dict(attr=attr, static=static, static_not_sequenced=True,
                      not_sequenced_val=self.ns_color.color if color else self.ns_number.value)
        if not color and values["not_sequenced_val"] is None:
            raise ValueError("Enter a value for not-sequenced observations.")

        if static:
            values["val"] = self.value_color.color if color else self.value_number.value
            if values["val"] is None:
                raise ValueError("Enter a value.")
            return values

        colname = self.colname_input.value.strip()
        if colname not in self._column_choices(attr):
            raise ValueError("Choose a categorical metadata column." if color
                             else "Choose a gene or a numeric metadata column.")
        values.update(colname=colname, colname_hovered=0 in self.hover_checkbox.active)
        if not color:
            low, high = self.min_input.value, self.max_input.value
            if low is None or high is None:
                raise ValueError("Enter the output range.")
            values.update(min_val=low, max_val=high)
        return values

    def _on_submit(self):
        try:
            values = self.retrieve_form_values()
        except ValueError as exc:
            self._error(str(exc))
            return
        editing = self._editing
        # add_parameter replaces the parameter already bound to the attribute
        self.component.add_parameter(**values)
        self.component.update_hover_tooltips()
        self._start_adding()
        done = "updated" if editing else "added"
        self.status.text = f"<p style='color:#5ee0c2;font-size:11px;margin:0'>{escape(values['attr'])} {done}.</p>"

    def _remove(self, attr):
        self.component.remove_parameter(attr)
        self._start_adding()
        self.status.text = f"<p style='color:#f0a500;font-size:11px;margin:0'>{escape(attr)} removed.</p>"
