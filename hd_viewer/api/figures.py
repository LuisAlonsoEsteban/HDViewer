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
                          CustomJSHover, HoverTool, LassoSelectTool, LinearColorMapper, TapTool, WMTSTileSource)
from bokeh.plotting import figure
from bokeh.themes import Theme
from bokeh.transform import transform

from ..board.layout.figure import FigurePanel
from ..core.parameters import is_color
from ..core.workspace import micrometres_per_unit
from ..startup.data.image.image import register_tiles, tile_extension


# Dark theme of the Bokeh UI, for figures shown in the panels (charts, QC, embeddings)
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


def map_figure(loaded, views, prefix="", view_box=None):
    """The tissue map; `view_box` (x0, x1, y0, y1) keeps the user's view when
    the figure is rebuilt (e.g. after a Parameters tab change)."""
    panel = FigurePanel()
    panel.create_figure(loaded.image_processor.max_dim)
    plot = panel.p
    if view_box is not None:
        x0, x1, y0, y1 = view_box
        plot.x_range.update(start=x0, end=x1)
        plot.y_range.update(start=y0, end=y1)
    plot.name = "map"
    plot.add_tools("box_select", "lasso_select")
    _image(plot, loaded.image_processor, loaded.image_processor.mem <= _tile_threshold(), prefix)

    bins = next((view.sample for view in views.values() if view.kind == "bins"), None)
    if bins is not None and loaded.transform is not None:
        um_per_px = micrometres_per_unit(bins)
        panel.measure_tool.set_transform(loaded.transform, um_per_px)
    # The shared figure starts with no active inspector: activate the layers' hover
    plot.toolbar.active_inspect = [_layer(plot, view) for view in views.values()]
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
    url = tile_template(processor, prefix)
    plot.image_url(url=[url.replace("{Z}", "0").replace("{X}", "0").replace("{Y}", "0")], x=0, y=size, w=size, h=size)
    plot.add_layout(BoxAnnotation(name="viewport_box", left=0, right=size, bottom=0, top=size,
                                  fill_alpha=0.1, fill_color="black", line_width=2, line_color="black"))
    return json_item(plot, "thumbnail")


def _tile_threshold():
    from .. import settings
    return settings.TILE_MEM_THRESHOLD
