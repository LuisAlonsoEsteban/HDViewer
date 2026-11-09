"""Bokeh figures of the React front end, serialised with `json_item`.

A json_item is a one-shot snapshot, so every model the browser will update is
declared here and *named*: per layer `<id>:source`, `<id>:mapper`, `<id>:filter`,
`<id>:renderer`, `<id>:hover`, `<id>:rois`; on the map `draft`, `roi_box`, `roi_lasso`,
`roi_tap`; plus `viewport_box` on the thumbnail. The browser
finds them with `doc.get_model_by_name` and changes only data and properties.
CustomJS hooks report ranges, selections and hovered rows to `window.hdv`.
"""
from pathlib import Path

from bokeh.embed import json_item
from bokeh.models import (BooleanFilter, BoxAnnotation, BoxSelectTool, CDSView, ColumnDataSource, CustomJS,
                          CustomJSHover, HoverTool, InlineStyleSheet, LassoSelectTool, LinearColorMapper, TapTool,
                          WMTSTileSource)
from bokeh.plotting import figure
from bokeh.themes import Theme
from bokeh.transform import transform

from .. import settings
from .measure import MeasureTool
from ..core.parameters import is_color
from ..core.workspace import micrometres_per_unit
from ..startup.data.image.image import register_tiles, tile_extension


# Dark theme for figures shown in the panels (charts, QC, embeddings)
THEME = Theme(filename=str(Path(__file__).resolve().parents[1] / "theme.yaml"))


def emit(event, payload="{}"):
    """CustomJS that reports `payload` to the browser bridge (window.hdv)."""
    return CustomJS(code=f"if (window.hdv) window.hdv.emit({event!r}, {payload})")


def tile_template(processor, prefix=""):
    """(URL template for {Z}/{X}/{Y}, extension) of a processed image's tiles."""
    key = register_tiles(processor.tiles_path)
    return f"{prefix}/api/tiles/{key}/{{Z}}/{{X}}/{{Y}}.{tile_extension(key)}"


def _image(plot, processor, single_image, prefix):
    url = tile_template(processor, prefix)
    if single_image:
        size = processor.max_dim
        plot.image_url(url=[url.replace("{Z}", "0").replace("{X}", "0").replace("{Y}", "0")],
                       x=0, y=size, w=size, h=size, name="image:url")
    else:
        renderer = plot.add_tile(WMTSTileSource(url=url, tile_size=processor.tile_size, max_zoom=processor.max_zoom,
                                     min_zoom=0, wrap_around=False, x_origin_offset=0, y_origin_offset=0,
                                     initial_resolution=processor.max_dim / processor.tile_size))
        renderer.name = "image:tiles"


_emit = emit


def _layer(plot, view):
    lid = view.id
    columns = {"x": [], "y": [], "alpha": [], "color": []}
    if view.kind == "spots":
        columns["radius"] = []
    source = ColumnDataSource(data=columns, name=f"{lid}:source")
    mapper = LinearColorMapper(palette=list(view.rows.get_param("color").palette), low=-0.5, high=0.5,
                               name=f"{lid}:mapper")
    view_ = CDSView(filter=BooleanFilter(name=f"{lid}:filter"))
    style = dict(source=source, view=view_, color=transform("color", mapper), alpha="alpha",
                 hover_color=view.describe()["hover_colour"], name=f"{lid}:renderer")
    # Parameters tab: each extra attribute reads its own column (colours through
    # their own mapper), overriding what `color`/`alpha` set for it
    for param in view.extra_parameters():
        attr = param.attr
        columns[attr] = []
        if is_color(attr):
            palette = list(param.palette)
            style[attr] = transform(attr, LinearColorMapper(palette=palette, low=-0.5, high=len(palette) - 0.5,
                                                            name=f"{lid}:mapper:{attr}"))
        else:
            style[attr] = dict(field=attr)
    source.data = dict(columns)
    constants = view.rows.glyph_constants
    if view.kind == "bins":
        sample = view.sample
        transform_ = getattr(sample, "spatial_transform", None)
        renderer = plot.rect("x", "y", width=constants.get("width", sample.native_width),
                             height=constants.get("height", sample.native_height),
                             angle=constants.get("angle", transform_.angle if transform_ else 0.0), **style)
    elif view.kind == "cells":
        renderer = plot.patches("x", "y", **style)
    else:
        renderer = plot.circle("x", "y", radius="radius", **style)
    renderer.visible = view.visible
    hover = HoverTool(
        name=f"{lid}:hover", renderers=[renderer],
        tooltips=[tuple(row) for row in view.tooltips()],
        formatters={"$index": CustomJSHover(code=f"return window.hdv ? window.hdv.hoverName({lid!r}, value) : ''")},
    )
    plot.add_tools(hover)
    source.selected.js_on_change("indices", _emit("selection", f"{{layer: {lid!r}, indices: cb_obj.indices}}"))
    # Saved ROIs of the layer (outlines in their colours)
    rois = ColumnDataSource(dict(xs=[], ys=[], color=[]), name=f"{lid}:rois")
    plot.multi_polygons("xs", "ys", source=rois, fill_alpha=.06, line_color="color", fill_color="color",
                        line_width=2, level="overlay", name=f"{lid}:roi_renderer")
    return hover


# Keeps the image's aspect ratio when the plot is resized or zoomed (one image
# pixel is square on screen)
_PRESERVE_ASPECT = """
    let scheduled = false
    let observer = null
    export default ({plot}) => {
        const schedule = () => {
            if (scheduled) return
            scheduled = true
            requestAnimationFrame(() => {
                try {
                    if (plot.document == null) { observer?.disconnect(); return }
                    const w = plot.inner_width, h = plot.inner_height
                    const xr = plot.x_range, yr = plot.y_range
                    const dx = xr.end-xr.start, dy = yr.end-yr.start
                    if (!(w > 0 && h > 0 && dx > 0 && dy > 0)) return
                    const units = Math.max(dx/w, dy/h)
                    const width = units*w, height = units*h
                    if (Math.abs(width-dx) <= dx*1e-7 && Math.abs(height-dy) <= dy*1e-7) return
                    const cx = (xr.start+xr.end)/2, cy = (yr.start+yr.end)/2
                    xr.setv({start:cx-width/2, end:cx+width/2})
                    yr.setv({start:cy-height/2, end:cy+height/2})
                } finally { scheduled = false }
            })
        }
        if (observer == null) {
            const view = Bokeh.index.find_one(plot)
            if (view != null) {
                observer = new ResizeObserver(schedule)
                observer.observe(view.el)
            }
        }
        schedule()
    }
"""


TOOLBAR_CSS = """
:host { background: rgba(20, 37, 31, 0.92); border-radius: 8px; padding: 4px 2px; }
.bk-tool-icon { filter: invert(88%) sepia(8%) saturate(300%) hue-rotate(95deg); }
.bk-tool-button:hover { background-color: rgba(94, 224, 194, 0.16); }
.bk-tool-button.bk-active { background-color: rgba(94, 224, 194, 0.22); }
.bk-divider { border-color: #2f5a49; }
"""


def base_figure(max_dim):
    """The empty map: image-sized ranges, pan/zoom tools and the ruler."""
    plot = figure(sizing_mode="stretch_both", min_height=settings.PLOT_MIN_HEIGHT,
                  x_range=(0, max_dim), y_range=(0, max_dim),
                  background_fill_color=settings.PLOT_BACKGROUND_COLOR, tools="pan,wheel_zoom,save,reset")
    plot.grid.visible = False
    plot.axis.visible = False
    plot.toolbar.logo = None
    plot.toolbar.autohide = False
    # One dark surface around the image (no white frame); light toolbar icons on it
    plot.border_fill_color = settings.PLOT_BACKGROUND_COLOR
    plot.outline_line_color = None
    plot.min_border = 0
    plot.toolbar.stylesheets = [InlineStyleSheet(css=TOOLBAR_CSS)]
    callback = CustomJS(args=dict(plot=plot), code=_PRESERVE_ASPECT)
    plot.js_on_event("document_ready", callback)
    for axis in (plot.x_range, plot.y_range):
        axis.js_on_change("start", callback)
        axis.js_on_change("end", callback)
    return plot, MeasureTool(plot)


def map_figure(loaded, views, prefix="", view_box=None):
    """The tissue map; `view_box` (x0, x1, y0, y1) keeps the user's view when
    the figure is rebuilt (e.g. after a Parameters tab change)."""
    plot, measure_tool = base_figure(loaded.image_processor.max_dim)
    if view_box is not None:
        x0, x1, y0, y1 = view_box
        plot.x_range.update(start=x0, end=x1)
        plot.y_range.update(start=y0, end=y1)
    plot.name = "map"
    plot.add_tools("box_select", "lasso_select")
    _image(plot, loaded.image_processor, loaded.image_processor.mem <= _tile_threshold(), prefix)

    bins = next((view.sample for view in views.values() if view.kind == "bins"), None)
    if bins is not None and loaded.transform is not None:
        measure_tool.set_transform(loaded.transform, micrometres_per_unit(bins))
    elif bins is not None and bins.native_width > 0:
        # No spatial transform: one native bin spans `bin_size` µm
        measure_tool.set_scale(bins.bin_size / bins.native_width, bins.bin_size / bins.native_height)
    for view in views.values():
        _layer(plot, view)
    # Hover tooltips start off (each one queries the server): the toolbar's
    # hover button turns them on, and the browser keeps that choice on rebuilds
    plot.toolbar.active_inspect = []
    # ROI being drawn (preview before saving). The drawing tools select nothing:
    # the browser arms one, reads the `geometry` / `tap` event and disarms it.
    draft = ColumnDataSource(dict(xs=[], ys=[]), name="draft")
    plot.multi_polygons("xs", "ys", source=draft, fill_alpha=.05, fill_color="#ffffff", line_color="#ffffff",
                        line_width=2, line_dash="dashed", level="overlay")
    plot.add_tools(BoxSelectTool(name="roi_box", renderers=[], continuous=False, description="Draw ROI bounding box"),
                   LassoSelectTool(name="roi_lasso", renderers=[], continuous=False, description="Draw ROI lasso"),
                   TapTool(name="roi_tap", renderers=[], description="Set ROI seed"))
    plot.js_on_event("selectiongeometry", _emit("geometry", "{geometry: cb_obj.geometry, final: cb_obj.final}"))
    plot.js_on_event("tap", _emit("tap", "{x: cb_obj.x, y: cb_obj.y}"))
    for axis in (plot.x_range, plot.y_range):
        for attr in ("start", "end"):
            axis.js_on_change(attr, _emit("ranges"))
    return json_item(plot, "map")


def thumbnail_figure(loaded, prefix=""):
    processor = loaded.image_processor
    size = processor.max_dim
    plot = figure(name="thumbnail", sizing_mode="scale_width", x_range=(0, size), y_range=(0, size),
                  match_aspect=True, toolbar_location=None, tools="")
    plot.grid.visible = False
    plot.axis.visible = False
    plot.background_fill_color = None
    plot.border_fill_color = None
    plot.outline_line_color = None
    plot.min_border = 0
    url = tile_template(processor, prefix)
    plot.image_url(url=[url.replace("{Z}", "0").replace("{X}", "0").replace("{Y}", "0")], x=0, y=size, w=size, h=size)
    plot.add_layout(BoxAnnotation(name="viewport_box", left=0, right=size, bottom=0, top=size,
                                  fill_alpha=0.12, fill_color="#5ee0c2", line_width=2, line_color="#5ee0c2"))
    return json_item(plot, "thumbnail")


def _tile_threshold():
    return settings.TILE_MEM_THRESHOLD
