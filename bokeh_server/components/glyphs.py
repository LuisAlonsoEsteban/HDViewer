
from bokeh.models import (
    HoverTool, 
    CheckboxGroup,
    AutocompleteInput,
    Slider,
    ColorPicker
)

from ..manager import SourceManager, DynamicSourceManager
from ..parameter import StaticParameter, Scaler

from ..loader import LazyLoader


class Component:

    def __init__(
            self,
            sample,
            lazy_loading=False,
            hover_color=None,
            alias='Spot',
            ):
        
        self.sample = sample
        self.lazy_loading = lazy_loading

        self.set_manager()

        self.hover_color = hover_color
        self.alias = alias

        self.glyph = None

    def set_manager(self):
        if self.lazy_loading:
            self.manager = DynamicSourceManager(self.sample)
        else:
            self.manager = SourceManager(self.sample)

    def attach_loader(self, p, doc):
        if self.lazy_loading and self.glyph is not None:
            loader = LazyLoader(self.manager, p)
            loader.add_range_callback(self, doc)

    def add_parameter(
            self, 
            attr='alpha', 
            static=False,
            val=0.5,
            min_val=0.0,
            max_val=1.0,
            colname=None,
            static_not_sequenced=True,
            not_sequenced_val=0.0,
            colname_hovered=False
            ):
        
        assert attr in ['color', 'alpha', 'fill_alpha', 'line_width', 'line_alpha'], "Incorrect attribute"
        
        if static:
            self.manager.add_parameter(
                StaticParameter(
                    val, 
                    attr=attr, 
                    not_sequenced_val=not_sequenced_val
                    )
            )
        else:
            self.manager.add_parameter(
                Scaler(
                    min_val,
                    max_val,
                    colname,
                    attr=attr,
                    static_not_sequenced=static_not_sequenced,
                    not_sequenced_val=not_sequenced_val,
                    colname_hovered=colname_hovered
                    )
            )


    def add_glyph(self, p):
        pass

    def get_hover_tooltips(self):
        items = ["x: @x", "y: @y"]
        for param in self.manager.parameters:
            if not param.is_static and param.colname_hovered:
                items.append(f"{param.colname}: @{param.colname}") # TODO: Add column data to source

        if len(items) > 0:
            return ", ".join(items)
        else:
            return ""
        
    def update_hover_tooltips(self):
        if hasattr(self, 'hover'):
            self.hover.tooltips = self.get_hover_tooltips()

    def add_hover(self, p):
        if self.hover_color is not None:
            self.hover = HoverTool(tooltips=[(self.alias, self.get_hover_tooltips())], renderers=[self.glyph])
            p.add_tools(self.hover)

    def set_visibility(self, visible:bool, disable_checkbox=False):
        self.glyph.visible = visible
        for widget in self.widgets:
            if widget is not self.checkbox or disable_checkbox:
                widget.disabled = not visible

    def get_visible_checkbox(self):
        self.checkbox = CheckboxGroup(labels=[f"Show {self.alias}s"], active=[0])

        def on_checkbox_change(attr, old, new):
            visible = 0 in new
            self.set_visibility(visible)

        self.checkbox.on_change('active', on_checkbox_change)

        return self.checkbox

    def get_colnames(self):
        return self.sample.get_colnames()
    
    def get_autocomplete_input(self, doc):
        autocomplete_input = AutocompleteInput(
                                    title="Gene/Metadata", 
                                    value=self.get_colnames()[0], 
                                    completions=self.get_colnames(), 
                                    min_characters=2, 
                                    case_sensitive=False, 
                                    restrict=True
                                    )
        
        def on_input_change(_, old, new):
            def update_colname():
                self.manager.update_colname(old, new)
                self.update_hover_tooltips()

            doc.add_next_tick_callback(update_colname)

        autocomplete_input.on_change('value', on_input_change)

        return autocomplete_input
    
    def get_widgets(self, doc):
        self.widgets = [self.get_visible_checkbox(), self.get_autocomplete_input(doc)]
        for param in self.manager.parameters:

            if param.attr == 'alpha' and not param.is_static:
                widget = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title=f"{self.alias}s alpha")

                def update_alpha(attr, old, new):
                    def set_new_max_alpha():
                        self.manager.update_param('alpha', new_maxv=new)
                    doc.add_next_tick_callback(set_new_max_alpha)

                widget.on_change('value', update_alpha)

            elif param.attr == 'color':
                widget = ColorPicker(title="Select Spots Color", color='#6baed6')

                def update_color(attr, old, new):
                    def set_new_color():
                        self.manager.update_param('color', new_val=new)
                    doc.add_next_tick_callback(set_new_color)

                widget.on_change('color', update_color)

            self.widgets.append(widget)

        return self.widgets
    

class Circles(Component):

    def __init__(self, radius, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius

    def add_glyph(self, p):
        self.glyph = p.circle(
                        x='x', 
                        y='y', 
                        radius=self.radius, 
                        hover_color=self.hover_color,
                        **{param.attr: param.attr for param in self.manager.parameters},
                        source=self.manager.source,
                    )
        
class Squares(Component):

    def __init__(self, width, height, angle, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.angle=angle

    def add_glyph(self, p):
        self.glyph = p.rect(
                        x='x', 
                        y='y', 
                        width=self.width, 
                        height=self.height,
                        angle=self.angle,
                        hover_color=self.hover_color,
                        **{param.attr: param.attr for param in self.manager.parameters},
                        source=self.manager.source,
                    )
        
class Patches(Component):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_glyph(self, p):
        self.glyph = p.patches(
                        xs='x', 
                        ys='y', 
                        hover_color=self.hover_color,
                        **{param.attr: param.attr for param in self.manager.parameters},
                        source=self.manager.source,
                    )



