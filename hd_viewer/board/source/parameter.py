"""Bokeh side of value encoding: the core parameters plus a LinearColorMapper
mirroring each colour palette (a glyph colour spec transform)."""
from bokeh.models import LinearColorMapper

from ...core.encoding import (            # noqa: F401  (re-exported)
    ColorCodes, DynamicParameter, Parameter, Scaler,
    CategoricalColor as _CategoricalColor, StaticParameter as _StaticParameter,
)


class MappedColors:
    """Keeps `mapper` in step with the core `palette`."""

    mapper = None

    def _palette_changed(self):
        palette = list(self.palette)
        if self.mapper is None:
            self.mapper = LinearColorMapper(palette=palette, low=-0.5, high=len(palette) - 0.5)
        elif list(self.mapper.palette) != palette:
            # Half-integer bounds: code k falls in the middle of palette bin k
            self.mapper.update(palette=palette, low=-0.5, high=len(palette) - 0.5)


class StaticParameter(MappedColors, _StaticParameter):
    pass


class CategoricalColor(MappedColors, _CategoricalColor):
    pass
