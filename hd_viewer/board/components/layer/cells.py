
from ._layer import Layer


class Patches(Layer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_glyph(self):
        self.glyph_renderer = self.figure.patches(
                                    xs='x', 
                                    ys='y',
                                    hover_color=self.hover_color,
                                    **self.glyph_specs(),
                                    source=self.manager.source,
                                    view=self.make_view(),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()