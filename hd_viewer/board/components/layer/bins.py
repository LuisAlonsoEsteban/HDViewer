from ._layer import Layer
from .utils.bin_info import BinInfoOverlay


class Squares(Layer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bin_info = BinInfoOverlay(self)

    def _constants(self):
        sample = self.sample
        transform = getattr(sample, "spatial_transform", None)
        defaults = dict(width=getattr(sample, "native_width", 1.0),
                        height=getattr(sample, "native_height", 1.0),
                        angle=transform.angle if transform else 0.0)
        return {**defaults, **getattr(self.manager, "glyph_constants", {})}

    def add_glyph(self):
        # Width, height and angle are the same for every bin of a display level:
        # glyph values rather than per-row columns
        self.glyph_renderer = self.figure.rect(
            x='x',
            y='y',
            **self._constants(),
            hover_color=self.hover_color,
            **self.glyph_specs(),
            source=self.manager.source,
            view=self.make_view(),
        )

        self.bin_info.bind_renderer(self.glyph_renderer)
        self.attach_selection_callbacks()
        self.apply_selected_obs_names()

    def apply_glyph_constants(self, constants):
        renderer = self.glyph_renderer
        if renderer is None or not constants:
            return
        for glyph in (renderer.glyph, renderer.selection_glyph, renderer.nonselection_glyph,
                      renderer.hover_glyph, renderer.muted_glyph):
            if glyph is None or isinstance(glyph, str):
                continue
            changes = {name: value for name, value in constants.items() if getattr(glyph, name, None) != value}
            if changes:
                glyph.update(**changes)
