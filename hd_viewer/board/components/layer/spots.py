
from ._layer import Layer


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
                                    **self.glyph_specs(),
                                    source=self.manager.source,
                                    view=self.make_view(),
                                )
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()