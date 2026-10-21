"""Zoom bar under the map, with the checkpoint from which cells are drawn."""

import math

from bokeh.models import Div

# Right end of the bar; the scale is logarithmic (1x = whole slide)
MAX_ZOOM = 128


class ZoomBar:
    """
    Shows how far the main plot is zoomed in, relative to the whole slide.

    With a cells gate (`set_cells_gate`), marks the zoom from which cells layers
    are loaded and drawn: cells cannot be aggregated like bins, so they are only
    shown once the view is small enough.
    """

    def __init__(self, figure, full_span):
        self.figure = figure
        self.full_span = float(full_span)   # x extent of the whole slide, figure units
        self._gate = None                   # callable -> max view width (figure units) for cells
        self._gate_label = ""
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

    def set_cells_gate(self, max_view_span, label=""):
        """max_view_span: callable giving the widest view (figure units) showing cells."""
        self._gate = max_view_span
        self._gate_label = label
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

    def gate_zoom(self):
        span = self._gate() if self._gate is not None else None
        return self.full_span / span if span else None

    def cells_visible(self):
        zoom, gate = self.zoom(), self.gate_zoom()
        return gate is None or (zoom is not None and zoom >= gate)

    @staticmethod
    def _position(zoom):
        """Percentage along the bar (log scale)."""
        return 100 * min(1.0, max(0.0, math.log(max(zoom, 1.0)) / math.log(MAX_ZOOM)))

    def render(self):
        zoom = self.zoom()
        if zoom is None:
            return
        gate = self.gate_zoom()
        fill = self._position(zoom)

        marker, status = "", ""
        if gate is not None:
            at = self._position(gate)
            marker = (f"<div title='Cells are drawn from {gate:.1f}× on' style='position:absolute;left:{at:.2f}%;"
                      f"top:-5px;bottom:-5px;width:2px;background:#f0c36b'></div>"
                      f"<div style='position:absolute;left:{at:.2f}%;top:-19px;transform:translateX(-50%);"
                      f"color:#f0c36b;font-size:10px;white-space:nowrap'>cells</div>")
            if zoom >= gate:
                status = "<span style='color:#5ee0c2'>Cells shown</span>"
            else:
                status = (f"<span style='color:#f0c36b'>Cells hidden · zoom in to {gate:.1f}×"
                          f"{' (' + self._gate_label + ')' if self._gate_label else ''}</span>")

        text = (
            "<div style='display:flex;align-items:center;gap:12px;font-size:12px;color:#b7cbd1'>"
            f"<span style='white-space:nowrap;min-width:78px'>Zoom <b style='color:#e5edf0'>{zoom:.1f}×</b></span>"
            "<div role='meter' aria-label='Zoom level' "
            f"aria-valuemin='1' aria-valuemax='{MAX_ZOOM}' aria-valuenow='{zoom:.1f}' "
            "style='position:relative;flex:1;height:8px;margin-top:12px;background:#394c52;border-radius:4px'>"
            f"<div style='width:{fill:.2f}%;height:100%;background:#5ee0c2;border-radius:4px'></div>"
            f"{marker}</div>"
            f"<span style='white-space:nowrap'>{status}</span>"
            "</div>"
        )
        if self.div.text != text:
            self.div.text = text
