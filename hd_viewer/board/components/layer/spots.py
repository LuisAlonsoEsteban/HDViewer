
from ._layer import Layer

from bokeh.models import CDSView, BooleanFilter

class Circles(Layer):

    def __init__(self, radius, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius

    def add_glyph(self):
        self.glyph_renderer = self.figure.circle(
                                    x='x', 
                                    y='y', 
                                    radius=self.radius, 
                                    hover_color=self.hover_color,
                                    **{param.attr: param.attr for param in self.manager.parameters},
                                    source=self.manager.source,
                                    view=CDSView(filter=BooleanFilter(booleans=self.view)),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()