"""Readout of the region currently visible on the main plot, in image coordinates."""

from html import escape
import numpy as np

from bokeh.layouts import row
from bokeh.models import Button, ColumnDataSource, CustomJS, Div, Node, Panel, SVGIcon

_COPY_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24"
fill="none" stroke="#e5edf0" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
<rect x="9" y="9" width="12" height="12" rx="2"/><path d="M5 15H4a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1h10a1 1 0 0 1 1 1v1"/></svg>"""

# Shared by the server (initial text) and the browser (live updates); keep in sync
_JS_FORMAT = """
function fmt_len(um) {
    return um >= 1000 ? `${(um / 1000).toFixed(2)} mm` : `${um.toFixed(0)} µm`
}
function roi(x0, x1, y0, y1, scale) {
    const px_x = scale.data.px_x[0], px_y = scale.data.px_y[0], um = scale.data.um_per_px[0]
    const scaled = px_x != null && px_y != null
    const fx = scaled ? px_x : 1, fy = scaled ? px_y : 1
    let box
    const inverse = scale.data.inverse?.[0]
    if (inverse != null) {
        const corners = [[x0,y0], [x1,y0], [x1,y1], [x0,y1]].map(([x,y]) =>
            [inverse[0]*x + inverse[1]*y + inverse[2], inverse[3]*x + inverse[4]*y + inverse[5]])
        box = [Math.min(...corners.map(p=>p[0])), Math.min(...corners.map(p=>p[1])),
               Math.max(...corners.map(p=>p[0])), Math.max(...corners.map(p=>p[1]))].map(Math.round)
    } else {
        box = [x0 * fx, y0 * fy, x1 * fx, y1 * fy].map(Math.round)
    }
    const n = (v) => v.toLocaleString("en-US")
    const lines = [`x  ${n(box[0])} – ${n(box[2])}`,
                   `y  ${n(box[1])} – ${n(box[3])}`]
    if (inverse != null) { lines.unshift("Original image · enclosing bounds") }
    if (scaled && um != null) {
        lines.push(`${fmt_len((box[2] - box[0]) * um)} × ${fmt_len((box[3] - box[1]) * um)}`)
    }
    return {lines, box}
}
"""


def _fmt_len(um):
    return f"{um / 1000:.2f} mm" if um >= 1000 else f"{um:.0f} µm"


class ROIReadout:
    """
    Shows the visible region of `figure` as a bounding box in full-resolution
    image pixels (Space Ranger's pxl_col/row_in_fullres; plot units until the
    scale is known) and its physical size. Updates live while panning/zooming;
    the copy button puts "x0,y0,x1,y1" on the clipboard. Call `place_on` to
    show it on a plot (e.g. the thumbnail).
    """

    STYLE = dict(font="11px/1.5 monospace", color="#e5edf0", white_space="pre", margin="0")

    def __init__(self, figure):
        self.figure = figure
        # Inverse affine calibration for original full-res image coordinates.
        # Per-axis scale fields retain support for scale-only callers.
        self.scale = ColumnDataSource(data=dict(px_x=[None], px_y=[None], um_per_px=[None], inverse=[None]))

        self.text = Div(styles=self.STYLE)
        self.copy_button = Button(
            label="", icon=SVGIcon(svg=_COPY_SVG, size=14), width=24, height=22, margin=0,
            align="start",
            stylesheets=["""
                .bk-btn { background: #23383e; border: 1px solid #5b8585;
                    border-radius: 6px; padding: 0; }
                .bk-btn:hover { background: #30464b; }
                .bk-btn:focus-visible { outline: 2px solid #5ee0c2; outline-offset: 2px; }
            """],
        )
        self.panel = None

        self._update = CustomJS(
            args=dict(fig=figure, scale=self.scale, div=self.text),
            code=_JS_FORMAT + """
                const {lines} = roi(fig.x_range.start, fig.x_range.end,
                                    fig.y_range.start, fig.y_range.end, scale)
                const text = lines.join("\\n")
                if (div.text != text) { div.text = text }
            """,
        )
        self.copy_button.js_on_click(CustomJS(
            args=dict(fig=figure, scale=self.scale),
            code=_JS_FORMAT + """
                const {box} = roi(fig.x_range.start, fig.x_range.end,
                                  fig.y_range.start, fig.y_range.end, scale)
                const text = box.join(",")
                // The clipboard API needs https or localhost; otherwise let the user copy
                if (navigator.clipboard && window.isSecureContext) {
                    navigator.clipboard.writeText(text)
                } else {
                    window.prompt("Visible region (x0,y0,x1,y1, full-res px):", text)
                }
            """,
        ))
        self.scale.js_on_change("change", self._update)
        # Ranges can be replaced (e.g. fitted to a low-res image): follow them
        figure.on_change("x_range", lambda attr, old, new: self._bind(new))
        figure.on_change("y_range", lambda attr, old, new: self._bind(new))
        self._bind(figure.x_range)
        self._bind(figure.y_range)

    def place_on(self, plot):
        """Show the readout in the top-right corner of `plot`."""
        self.panel = Panel(
            name="roi_readout",
            position=Node(target="frame", symbol="top_right"),
            anchor="top_right",
            elements=[row(self.text, self.copy_button, margin=0, spacing=8)],
            styles=dict(background="rgba(24,39,45,0.94)", border="1px solid #5b8585",
                        border_radius="8px", padding="6px 8px", margin_left="-8px",
                        margin_top="8px", pointer_events="auto"),
        )
        plot.elements.append(self.panel)

    def _bind(self, rng):
        rng.js_on_change("start", self._update)
        rng.js_on_change("end", self._update)
        self.render()

    def set_scale(self, px_per_unit_x, px_per_unit_y, um_per_px=None):
        self.scale.data = dict(px_x=[float(px_per_unit_x)], px_y=[float(px_per_unit_y)],
                               um_per_px=[None if um_per_px is None else float(um_per_px)], inverse=[None])
        self.render()

    def set_transform(self, transform, um_per_px=None):
        self.scale.data = dict(px_x=[1.], px_y=[1.], um_per_px=[um_per_px],
                               inverse=[transform.inverse.ravel().tolist()])
        self.render()

    def bbox(self):
        """Visible region as (x0, y0, x1, y1), in full-res pixels once scaled."""
        x_range, y_range = self.figure.x_range, self.figure.y_range
        inverse = self.scale.data.get('inverse', [None])[0]
        if inverse is not None:
            matrix = np.asarray(inverse).reshape(3, 3)
            corners = np.array([[x_range.start, y_range.start], [x_range.end, y_range.start],
                                [x_range.end, y_range.end], [x_range.start, y_range.end]])
            points = corners @ matrix[:2, :2].T + matrix[:2, 2]
            return tuple(round(v) for v in (*points.min(axis=0), *points.max(axis=0)))
        fx = self.scale.data["px_x"][0] or 1
        fy = self.scale.data["px_y"][0] or 1
        return tuple(round(v) for v in (x_range.start * fx, y_range.start * fy,
                                        x_range.end * fx, y_range.end * fy))

    def render(self):
        """Server-side text, so the readout is right before the first pan."""
        if self.figure.x_range.start is None or self.figure.y_range.start is None:
            return
        scaled = self.scale.data["px_x"][0] is not None
        x0, y0, x1, y1 = self.bbox()
        lines = [f"x  {x0:,} – {x1:,}", f"y  {y0:,} – {y1:,}"]
        if self.scale.data.get("inverse", [None])[0] is not None:
            lines.insert(0, "Original image · enclosing bounds")
        um = self.scale.data["um_per_px"][0]
        if scaled and um is not None:
            lines.append(f"{_fmt_len((x1 - x0) * um)} × {_fmt_len((y1 - y0) * um)}")
        self.text.text = escape("\n".join(lines))
