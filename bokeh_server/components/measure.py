"""Ruler tool: click two points on the plot to measure the distance between them."""

import base64

from bokeh.models import (
    ColumnDataSource,
    CustomAction,
    CustomJS,
    Label,
    PointDrawTool,
)

_RULER_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
stroke="#555" stroke-width="1.6" stroke-linecap="round">
<rect x="2" y="8" width="20" height="8" rx="1" transform="rotate(-30 12 12)"/>
<path d="M7.2 14.2l-1-1.7M10 12.6l-1.5-2.6M12.8 11l-1-1.7M15.6 9.4l-1.5-2.6M18.4 7.8l-1-1.7"/></svg>"""

_CLEAR_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none"
stroke="#555" stroke-width="1.8" stroke-linecap="round">
<path d="M4 20h16M7 16l9-9M8 8l8 8"/></svg>"""


def _icon(svg):
    return "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()


class MeasureTool:
    """
    Two-point ruler for the main figure.

    Activate it from the toolbar, click the start and end points, and drag
    either end to adjust; a third click starts a new measurement. The distance
    is shown in µm once `set_scale` has been called (figure units otherwise).
    """

    COLOR = "#ffd166"

    def __init__(self, figure):
        self.figure = figure

        self.points = ColumnDataSource(data=dict(x=[], y=[]))
        # µm per figure unit along each axis: tiles stretch the image to a square,
        # so the two axes can have different scales
        self.scale = ColumnDataSource(data=dict(um_x=[None], um_y=[None]))

        line = figure.line('x', 'y', source=self.points, line_color=self.COLOR,
                           line_width=2, line_dash="dashed")
        ends = figure.scatter('x', 'y', source=self.points, size=10, marker="circle",
                              fill_color=self.COLOR, line_color="#1b1b1b")
        line.level = ends.level = "overlay"

        self.label = Label(
            x=0, y=0, text="", visible=False, level="overlay",
            x_offset=8, y_offset=8,
            text_color="#1b1b1b", text_font_size="13px", text_font_style="bold",
            background_fill_color=self.COLOR, background_fill_alpha=0.9,
            border_line_color="#1b1b1b", border_radius=4, padding=4,
        )
        figure.add_layout(self.label)

        self.draw_tool = PointDrawTool(
            renderers=[ends], num_objects=2, add=True, drag=True,
            description="Measure distance (click two points, drag to adjust)",
            icon=_icon(_RULER_SVG),
        )
        clear_tool = CustomAction(
            description="Clear measurement",
            icon=_icon(_CLEAR_SVG),
            callback=CustomJS(args=dict(points=self.points),
                              code="points.data = {x: [], y: []}"),
        )
        figure.add_tools(self.draw_tool, clear_tool)

        # 'change' (not 'data') also fires while a point is being dragged
        update = CustomJS(args=dict(points=self.points, scale=self.scale, label=self.label), code="""
            const xs = points.data.x, ys = points.data.y
            if (xs.length < 2) { label.visible = false; return }
            const [x0, x1] = [xs[xs.length - 2], xs[xs.length - 1]]
            const [y0, y1] = [ys[ys.length - 2], ys[ys.length - 1]]
            const um_x = scale.data.um_x[0], um_y = scale.data.um_y[0]
            let text
            if (um_x != null && um_y != null) {
                const d = Math.hypot((x1 - x0) * um_x, (y1 - y0) * um_y)
                text = d >= 1000 ? `${(d / 1000).toFixed(2)} mm` : `${d.toFixed(d < 10 ? 2 : 1)} µm`
            } else {
                text = `${Math.hypot(x1 - x0, y1 - y0).toFixed(1)} px`
            }
            label.setv({x: (x0 + x1) / 2, y: (y0 + y1) / 2, text: text, visible: true})
        """)
        self.points.js_on_change('change', update)
        self.scale.js_on_change('change', update)

    def set_scale(self, um_per_unit_x, um_per_unit_y):
        self.scale.data = dict(um_x=[float(um_per_unit_x)], um_y=[float(um_per_unit_y)])
