
from bokeh.models import (
    Dropdown,
    CheckboxGroup,
    AutocompleteInput,
    TextInput,
    Button
)
from bokeh.layouts import column

class ParamForm:

    def __init__(self, glyph):
        self.glyph = glyph

    def get_layout(self):
        self.dropdown = Dropdown(label="Attribute", menu=self.glyph.get_available_attrs())
        self.static_checkbox = CheckboxGroup(labels=["Static"], active=[])
        self.colname_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=None, 
                                    completions=self.glyph.get_colnames(), 
                                    min_characters=2, 
                                    case_sensitive=False, 
                                    restrict=True
                                    )
        self.val_input = TextInput(value="", title="Value")
        self.min_val_input = TextInput(value="", title="Max value")
        self.max_val_input = TextInput(value="", title="Min value")
        self.sns_checkbox = CheckboxGroup(labels=["Static Not-Sequenced"], active=[])
        self.not_sequenced_val_input = TextInput(value="", title="Not-sequenced value")
        self.colname_hovered_checkbox = CheckboxGroup(labels=["Static Not-Sequenced"], active=[])

        self.confirm_button = Button(label="Confirm")

        def on_click_button():
            param_attrs = {
                'attr': self.dropdown.value,
                'is_static': 0 in self.static_checkbox.active,
                'colname': self.colname_input.value,
                'val': self.val_input.value,
                'min_val': self.min_val_input.value,
                'max_val': self.max_val_input.value,
                'static_not_sequenced': 0 in self.sns_checkbox.active,
                'not_sequenced_val': self.not_sequenced_val_input.value,
                'colname_hovered': 0 in self.colname_hovered_checkbox.active 
            }
            print(param_attrs)
            # self.glyph.add_parameter(**param_attrs)

        self.confirm_button.on_click(on_click_button)

        return column(
            self.dropdown,
            self.static_checkbox,
            self.colname_input,
            self.val_input,
            self.min_val_input,
            self.max_val_input,
            self.sns_checkbox,
            self.not_sequenced_val_input,
            self.colname_hovered_checkbox,
            self.confirm_button
        )