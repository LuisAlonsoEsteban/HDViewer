"""The loaded workspace: figures, and per-layer data as binary frames.

Viewport, style and state responses are frames (api.frames) holding the
operations to apply to the layer's BokehJS data source, with `base`/`version`
data versions: a browser whose version differs from `base` asks for /state.
"""
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import FileResponse, Response
from starlette.background import BackgroundTask
from pydantic import BaseModel

from ... import settings
from ...core.export import export_layer
from ...core.workspace import micrometres_per_unit
from .. import figures, frames
from ..views import build_views
from . import current_session

router = APIRouter()
FRAME = "application/octet-stream"


def open_workspace(session):
    loaded = session.loaded
    bins = next((sample for kind, sample, _ in loaded.samples.values() if kind == "bins"), None)
    session.views = build_views(loaded, micrometres_per_unit(bins) if bins is not None else None)


def _loaded(session):
    if session.loaded is None:
        raise HTTPException(409, "No dataset is loaded.")
    return session.loaded


def _view(session, layer_id):
    _loaded(session)
    view = session.views.get(layer_id)
    if view is None:
        raise HTTPException(404, f"Unknown layer: {layer_id}")
    return view


def _frame(view, ops):
    return Response(frames.encode(*view.frame_parts(ops)), media_type=FRAME)


def _prefix(request):
    return request.scope.get("root_path", "")


def _inverse(transform, bins):
    """Map units -> full-resolution image pixels (3x3, row-major): the spatial
    transform's inverse, or (no transform) the image's tile scale per axis."""
    if transform is not None:
        return transform.inverse.ravel().tolist()
    if bins is None or not getattr(bins, "tile_scale_xfactor", None) or not getattr(bins, "tile_scale_yfactor", None):
        return None
    return [1 / bins.tile_scale_xfactor, 0.0, 0.0, 0.0, 1 / bins.tile_scale_yfactor, 0.0, 0.0, 0.0, 1.0]


@router.get("/workspace")
def workspace(request: Request, session=Depends(current_session)):
    loaded = _loaded(session)
    processor = loaded.image_processor
    bins = next((view.sample for view in session.views.values() if view.kind == "bins"), None)
    transform = loaded.transform
    return dict(
        layers=[view.describe() for view in session.views.values()],
        image=dict(size=loaded.image_processor.max_dim,
                   inverse=_inverse(transform, bins),
                   um_per_px=None if bins is None else micrometres_per_unit(bins)),
        tiles=dict(template=figures.tile_template(processor, _prefix(request)),
                   single=processor.mem <= settings.TILE_MEM_THRESHOLD, max_zoom=processor.max_zoom,
                   tile_size=processor.tile_size),
        reference=loaded.reference_name,
    )


@router.get("/figures/map")
def map_figure(request: Request, x0: float | None = None, x1: float | None = None, y0: float | None = None,
               y1: float | None = None, session=Depends(current_session)):
    """The map figure; with x0..y1 it opens at that view (rebuilds keep the view)."""
    box = (x0, x1, y0, y1)
    view_box = box if None not in box and x1 > x0 and y1 > y0 else None
    return figures.map_figure(_loaded(session), session.views, _prefix(request), view_box)


@router.get("/figures/thumbnail")
def thumbnail_figure(request: Request, session=Depends(current_session)):
    return figures.thumbnail_figure(_loaded(session), _prefix(request))


class ViewportRequest(BaseModel):
    seq: int
    x0: float
    x1: float
    y0: float
    y1: float
    width: float = 800
    height: float = 800


@router.post("/layers/{layer_id}/viewport")
def viewport(layer_id: str, body: ViewportRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        if body.seq <= view.seq:
            raise HTTPException(409, "stale request")
        view.seq = body.seq
        ops = view.update_viewport(body.x0, body.x1, body.y0, body.y1,
                                   max(body.width, 1.0), max(body.height, 1.0))
        return _frame(view, ops)


@router.get("/layers/{layer_id}/state")
def state(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return Response(frames.encode(*view.full_state()), media_type=FRAME)


class StyleRequest(BaseModel):
    visible: bool | None = None
    colour_by: str | None = None          # gene or metadata column
    scale: str | None = None              # counts | normalized | log1p
    colour: str | None = None             # hex colour for continuous values
    alpha: float | None = None            # maximum opacity
    aggregation: int | None = None        # bins: fixed factor; 0 follows the zoom


@router.patch("/layers/{layer_id}/style")
def style(layer_id: str, body: StyleRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    ops = []
    with view.lock:
        try:
            if body.visible is not None:
                ops += view.set_visible(body.visible)
            if body.scale is not None:
                ops += view.set_scale(body.scale)
            if body.colour_by is not None:
                ops += view.colour_by_column(body.colour_by)
            if body.colour is not None:
                ops += view.set_colour(body.colour)
            if body.alpha is not None:
                ops += view.set_alpha(min(max(body.alpha, 0.0), 1.0))
            if body.aggregation is not None:
                ops += view.set_aggregation(body.aggregation or None)
        except KeyError as exc:
            raise HTTPException(422, f"Unknown gene or metadata column: {exc.args[0]}")
        except ValueError as exc:
            raise HTTPException(422, str(exc))
        return _frame(view, ops)


@router.get("/layers/{layer_id}/columns")
def columns(layer_id: str, q: str = "", limit: int = Query(50, le=500), session=Depends(current_session)):
    return _view(session, layer_id).column_names(q, limit)


class RowsRequest(BaseModel):
    row_version: int
    rows: list[int]


@router.post("/layers/{layer_id}/names")
def names(layer_id: str, body: RowsRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return dict(row_version=view.rows.row_version, rows=view.names(body.row_version, body.rows[:1000]))


@router.post("/layers/{layer_id}/selection")
def select(layer_id: str, body: RowsRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        try:
            return dict(selected=view.select(body.row_version, body.rows))
        except ValueError:
            raise HTTPException(409, "stale rows")


@router.get("/layers/{layer_id}/selection")
def selection(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return dict(row_version=view.rows.row_version, rows=view.selected_positions(),
                    selected=len(view.selected))


# -- Parameters tab: extra glyph attributes ----------------------------------------------

class ParameterRequest(BaseModel):
    mode: str                                 # static | data
    value: str | float | None = None          # static: a hex colour, or a number
    column: str | None = None                 # data: a column of the layer
    min: float | None = None                  # data, numeric attributes: output range
    max: float | None = None
    not_sequenced: str | float | None = None  # value of not-sequenced rows
    hovered: bool = False                     # data: show the column on hover


@router.get("/layers/{layer_id}/parameters")
def parameters(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return view.describe_parameters()


@router.get("/layers/{layer_id}/parameters/{attr}/columns")
def parameter_columns(layer_id: str, attr: str, q: str = "", limit: int = Query(50, le=500),
                      session=Depends(current_session)):
    from ...core.parameters import EXTRA_ATTRS, column_choices
    view = _view(session, layer_id)
    if attr not in EXTRA_ATTRS:
        raise HTTPException(404, f"Unknown attribute: {attr}")
    query = q.lower()
    names = column_choices(view.sample, attr)
    return [name for name in names if not query or name.lower().startswith(query)][:limit]


@router.put("/layers/{layer_id}/parameters/{attr}")
def set_parameter(layer_id: str, attr: str, body: ParameterRequest, session=Depends(current_session)):
    """Changes the glyph: the browser then reloads the map figure (GET /figures/map)."""
    view = _view(session, layer_id)
    with view.lock:
        try:
            return view.set_parameter(attr, mode=body.mode, value=body.value, column=body.column,
                                      minv=body.min, maxv=body.max, not_sequenced=body.not_sequenced,
                                      hovered=body.hovered)
        except ValueError as exc:
            raise HTTPException(422, str(exc))


@router.delete("/layers/{layer_id}/parameters/{attr}", status_code=204)
def remove_parameter(layer_id: str, attr: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        try:
            view.remove_parameter(attr)
        except ValueError as exc:
            raise HTTPException(404, str(exc))
    return Response(status_code=204)


# -- downloads -------------------------------------------------------------------------------

@router.get("/layers/{layer_id}/export")
def export(layer_id: str, session=Depends(current_session)):
    """The layer's AnnData (.h5ad) with its clusters, groups and deconvolution
    results, written to a private temporary file removed once sent."""
    view = _view(session, layer_id)
    try:
        path = export_layer(view.sample, session.directory / "exports")
    except Exception as exc:
        raise HTTPException(500, f"The data could not be exported: {exc}")
    return FileResponse(path, media_type="application/x-hdf5", filename=f"{layer_id}.h5ad",
                        background=BackgroundTask(path.unlink, missing_ok=True))
