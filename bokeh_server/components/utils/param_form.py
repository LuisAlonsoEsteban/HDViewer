
from bokeh.models import (
    Dropdown,
    CheckboxGroup,
    AutocompleteInput,
    Button
)

class ParamForm:

    def __init__(self, available_attrs, colname_pool):
        self.available_attrs = available_attrs
        self.colname_pool = colname_pool

    def get_layout(self):
        self.dropdown = Dropdown(label="Attribute", menu=self.available_attrs)
        self.static_checkbox = CheckboxGroup(labels=["Static"], active=[])
        self.colname_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=None, 
                                    completions=self.colname_pool, 
                                    min_characters=2, 
                                    case_sensitive=False, 
                                    restrict=True
                                    )
        self.confirm_button = Button(label="Confirm")

        def on_click_button():
            return {
                'attr': self.dropdown.value,
                'is_static': 0 in self.static_checkbox.active,
                'colname': self.colname_input.value
                'val': ,
                'min_val': ,
                'max_val': ,
                'static_not_sequenced': ,
                'not_sequenced_val': ,
                'colname_hovered': 
            }