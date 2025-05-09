
from bokeh.plotting import column
from bokeh.models import (
    HoverTool, 
    AutocompleteInput,
    Slider,
    ColorPicker
)

from .manager import SourceManager, DynamicSourceManager
from .parameter import StaticParameter, Scaler

class Component:

    def __init__(
            self,
            sample,
            lazy_loading=False,
            hover_color=None,
            alias='Spot',
            visible=True
            ):
        
        self.sample = sample
        self.lazy_loading = lazy_loading

        self.set_manager()
        self.parameters = []

        self.hover_color = hover_color
        self.alias = alias

        self.glyph = None
        self.visible = visible

    def set_manager(self):
        if self.lazy_loading:
            self.manager = DynamicSourceManager(self.sample)
        else:
            self.manager = SourceManager(self.sample)
        self.parameters = self.manager.parameters

    def add_parameter(
            self, 
            attr='alpha', 
            static=False,
            val=0.5,
            min_val=0.0,
            max_val=1.0,
            colname=None,
            static_not_sequenced=None,
            not_sequenced_val=0.0
            ):
        
        assert attr in ['color', 'alpha', 'fill_alpha', 'line_width', 'line_alpha'], "Incorrect attribute"
        
        if static:
            self.manager.parameters.append(
                StaticParameter(
                    val, 
                    attr=attr, 
                    not_sequenced_val=not_sequenced_val
                    )
            )
        else:
            self.manager.parameters.append(
                Scaler(
                    min_val,
                    max_val,
                    colname,
                    attr=attr,
                    static_not_sequenced=static_not_sequenced,
                    not_sequenced_val=not_sequenced_val
                    )
            )


    def add_glyph(self, p):
        pass

    def get_hover_tooltips(self):
        items = []
        for param in self.parameters:
            if not param.is_static and param.colname_hovered:
                items.append(f"{param.colname}: {param.colname}")

        if len(items) > 0:
            return ", ".join(items)
        else:
            return ""

    def add_hover(self, p):
        if self.hover_color is not None:
            hover = HoverTool(tooltips=[(self.alias, self.get_hover_tooltips())], renderers=[self.glyph])
            p.add_tools(hover)

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
                if self.lazy_loading:
                    pass
                else:
                    self.manager.fetch_data()

            doc.add_next_tick_callback(update_colname)

        autocomplete_input.on_change('value', on_input_change)

        return autocomplete_input
    
    def get_widgets(self, doc):
        widgets = [self.get_autocomplete_input(doc)]
        for param in self.parameters:

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

            widgets.append(widget)

        return widgets

        
        
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
                        **{param.attr: param.attr for param in self.parameters},
                        source=self.source,
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
                        **{param.attr: param.attr for param in self.parameters},
                        source=self.source,
                    )
        
class Patches(Component):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_glyph(self, p):
        self.glyph = p.patches(
                        x='x', 
                        y='y', 
                        hover_color=self.hover_color,
                        **{param.attr: param.attr for param in self.parameters},
                        source=self.source,
                    )



