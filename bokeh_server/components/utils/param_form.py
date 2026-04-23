
from bokeh.models import (
    Select, 
    CheckboxGroup, 
    AutocompleteInput, 
    TextInput, 
    Button
)
from bokeh.layouts import column

class ParamForm:

    def __init__(self, component):
        self.component = component

    def get_layout(self, doc):
        self.attr_select = Select(title="Attribute", options=self.component.get_available_attrs())
        self.static_checkbox = CheckboxGroup(labels=["Static"], active=[])
        self.colname_input = AutocompleteInput(
            title="Gene/Metadata",
            value="",
            completions=self.component.get_colnames(),
            min_characters=2,
            case_sensitive=False,
            restrict=True
        )
        self.val_input = TextInput(value="", title="Value")
        self.min_val_input = TextInput(value="", title="Min Value")
        self.max_val_input = TextInput(value="", title="Max Value")
        self.sns_checkbox = CheckboxGroup(labels=["Static Not-Sequenced"], active=[])
        self.not_sequenced_val_input = TextInput(value="", title="Not-Sequenced Value")
        self.colname_hovered_checkbox = CheckboxGroup(labels=["Column Name Hovered"], active=[])

        self.confirm_button = Button(label="Confirm")

        # Setup JavaScript to handle dynamic enabling/disabling
        self.setup_callbacks(doc)

        return column(
            self.attr_select,
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
    
    def retrieve_form_values(self):
        attr = self.attr_select.value
        if 'color' in attr:
            dtype = lambda x: x
        else:
            dtype = lambda x: float(x) if x != '' else 0.0

        return {
                'attr': attr,
                'static': 0 in self.static_checkbox.active,
                'val': dtype(self.val_input.value),
                'min_val': dtype(self.min_val_input.value),
                'max_val': dtype(self.max_val_input.value),
                'colname': self.colname_input.value,
                'static_not_sequenced': 0 in self.sns_checkbox.active,
                'not_sequenced_val': dtype(self.not_sequenced_val_input.value),
                'colname_hovered': 0 in self.colname_hovered_checkbox.active
            }

    def setup_callbacks(self, doc):
        
        def on_click_button():
            param_attrs = self.retrieve_form_values()

            def update_parameter():
                self.component.add_parameter(**param_attrs)
                self.component.refresh_layout(doc)

            doc.add_next_tick_callback(update_parameter)

        self.confirm_button.on_click(on_click_button)
