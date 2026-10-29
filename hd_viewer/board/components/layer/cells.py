
from ._layer import Layer

from bokeh.models import CDSView, BooleanFilter

class Patches(Layer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_glyph(self):
        self.glyph_renderer = self.figure.patches(
                                    xs='x', 
                                    ys='y',
                                    hover_color=self.hover_color,
                                    **{param.attr: param.attr for param in self.manager.parameters},
                                    source=self.manager.source,
                                    view=CDSView(filter=BooleanFilter(booleans=self.view)),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()