
from ._layer import Layer
from .utils.bin_info import BinInfoOverlay

from bokeh.models import CDSView, BooleanFilter

class Squares(Layer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bin_info = BinInfoOverlay(self)

    def add_glyph(self):
        self.glyph_renderer = self.figure.rect(
            x='x',
            y='y',
            width='width',
            height='height',
            angle='angle',
            hover_color=self.hover_color,
            **{param.attr: param.attr for param in self.manager.parameters},
            source=self.manager.source,
            view=CDSView(filter=BooleanFilter(booleans=self.view)),
        )

        self.bin_info.bind_renderer(self.glyph_renderer)
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()