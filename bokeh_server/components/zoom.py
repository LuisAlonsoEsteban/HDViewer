"""Zoom bar under the map."""

import math

from bokeh.models import Div

# Right end of the bar; the scale is logarithmic (1x = whole slide)
MAX_ZOOM = 128


class ZoomBar:
    """
    Shows how far the main plot is zoomed in, relative to the whole slide.
    """

    def __init__(self, figure, full_span):
        self.figure = figure
        self.full_span = float(full_span)   # x extent of the whole slide, figure units
        self.div = Div(sizing_mode="stretch_width", margin=(6, 10, 4, 10))
        # Ranges can be replaced (e.g. fitted to a low-res image): follow them
        figure.on_change("x_range", lambda attr, old, new: self._bind(new))
        self._bind(figure.x_range)

    def _bind(self, rng):
        rng.on_change("start", lambda attr, old, new: self.render())
        rng.on_change("end", lambda attr, old, new: self.render())
        self.render()

    def set_full_span(self, span):
        self.full_span = float(span)
        self.render()

    # ------------------------------------------------------------------

    def view_span(self):
        x_range = self.figure.x_range
        if x_range.start is None or x_range.end is None:
            return None
        return x_range.end - x_range.start

    def zoom(self):
        span = self.view_span()
        return self.full_span / span if span and span > 0 else None

    @staticmethod
    def _position(zoom):
        """Percentage along the bar (log scale)."""
        return 100 * min(1.0, max(0.0, math.log(max(zoom, 1.0)) / math.log(MAX_ZOOM)))

    def render(self):
        zoom = self.zoom()
        if zoom is None:
            return
        fill = self._position(zoom)

        text = (
            "<div style='display:flex;align-items:center;gap:12px;font-size:12px;color:#b7cbd1'>"
            f"<span style='white-space:nowrap;min-width:78px'>Zoom <b style='color:#e5edf0'>{zoom:.1f}×</b></span>"
            "<div role='meter' aria-label='Zoom level' "
            f"aria-valuemin='1' aria-valuemax='{MAX_ZOOM}' aria-valuenow='{zoom:.1f}' "
            "style='position:relative;flex:1;height:8px;background:#394c52;border-radius:4px'>"
            f"<div style='width:{fill:.2f}%;height:100%;background:#5ee0c2;border-radius:4px'></div>"
            "</div></div>"
        )
        if self.div.text != text:
            self.div.text = text
