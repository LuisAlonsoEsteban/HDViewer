
from . import Panel

class ThumbnailPanel(Panel):

    thumb_p = None

    def __init__(self):
        super().__init__(name="thumbnail_plot")

    def populate(self, image_component, p, roi_readout=None):
        self.thumb_p = image_component.generate_thumbnail()

        image_component.add_lowres_image(self.thumb_p)
        image_component.add_viewport(self.thumb_p)
        image_component.add_viewport_callback(p)
        
        if roi_readout is not None and roi_readout.panel is None:
            roi_readout.place_on(self.thumb_p)

        self.set_children([self.thumb_p])


import numpy as np

from bokeh.models import Button, ColumnDataSource, CustomJS, Node, Panel, SVGIcon

from ..utils.icons import _COPY_SVG

# Visible region of the map as (x0, y0, x1, y1) in full-resolution image pixels;
# mirrors ViewportReadout.bbox.
_JS_BBOX = """
function viewport_box(x0, x1, y0, y1, scale) {
    const inverse = scale.data.inverse?.[0]
    if (inverse != null) {
        const corners = [[x0,y0], [x1,y0], [x1,y1], [x0,y1]].map(([x,y]) =>
            [inverse[0]*x + inverse[1]*y + inverse[2], inverse[3]*x + inverse[4]*y + inverse[5]])
        return [Math.min(...corners.map(p=>p[0])), Math.min(...corners.map(p=>p[1])),
                Math.max(...corners.map(p=>p[0])), Math.max(...corners.map(p=>p[1]))].map(Math.round)
    }
    const fx = scale.data.px_x[0] ?? 1, fy = scale.data.px_y[0] ?? 1
    return [x0 * fx, y0 * fy, x1 * fx, y1 * fy].map(Math.round)
}
"""


class ViewportReadout:
    """
    A button that copies the visible region of `figure` as `x0,y0,x1,y1` in
    full-resolution image pixels (Space Ranger's pxl_col/row_in_fullres; plot
    units until the scale is known). Nothing is displayed or recomputed while
    panning: the bounds are read only when the button is clicked.
    Call `place_on` to show it on a plot (e.g. the thumbnail).
    """

    LABEL = "Copy viewport coordinates"

    def __init__(self, figure):
        self.figure = figure
        # Inverse affine calibration for original full-res image coordinates.
        # Per-axis scale fields retain support for scale-only callers.
        self.scale = ColumnDataSource(data=dict(px_x=[None], px_y=[None], um_per_px=[None], inverse=[None]))

        self.copy_button = Button(
            label=self.LABEL, icon=SVGIcon(svg=_COPY_SVG, size=14), height=24, margin=0,
            stylesheets=["""
                .bk-btn { background: #23383e; border: 1px solid #5b8585; color: #e5edf0;
                    border-radius: 6px; padding: 2px 8px; font-size: 11px; }
                .bk-btn:hover { background: #30464b; }
                .bk-btn:focus-visible { outline: 2px solid #5ee0c2; outline-offset: 2px; }
            """],
        )
        self.panel = None

        self.copy_button.js_on_click(CustomJS(
            args=dict(fig=figure, scale=self.scale, button=self.copy_button, label=self.LABEL),
            code=_JS_BBOX + """
                const text = viewport_box(fig.x_range.start, fig.x_range.end,
                                          fig.y_range.start, fig.y_range.end, scale).join(",")
                // The clipboard API needs https or localhost; otherwise let the user copy
                if (navigator.clipboard && window.isSecureContext) {
                    navigator.clipboard.writeText(text).then(() => {
                        button.label = "Copied"
                        setTimeout(() => { button.label = label }, 1200)
                    })
                } else {
                    window.prompt("Visible region (x0,y0,x1,y1, full-res px):", text)
                }
            """,
        ))

    def place_on(self, plot):
        """Show the button in the top-right corner of `plot`."""
        self.panel = Panel(
            name="roi_readout",
            position=Node(target="frame", symbol="top_right"),
            anchor="top_right",
            elements=[self.copy_button],
            styles=dict(margin_left="-8px", margin_top="8px", pointer_events="auto"),
        )
        plot.elements.append(self.panel)

    def set_scale(self, px_per_unit_x, px_per_unit_y, um_per_px=None):
        self.scale.data = dict(px_x=[float(px_per_unit_x)], px_y=[float(px_per_unit_y)],
                               um_per_px=[None if um_per_px is None else float(um_per_px)], inverse=[None])

    def set_transform(self, transform, um_per_px=None):
        self.scale.data = dict(px_x=[1.], px_y=[1.], um_per_px=[um_per_px],
                               inverse=[transform.inverse.ravel().tolist()])

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
