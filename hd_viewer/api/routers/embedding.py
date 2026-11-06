"""2-D embeddings (PCA, UMAP, t-SNE) of a layer's native observations.

An embedding is a background job (events `job` with name "embedding"; it shares
the layer's slot with clustering and DE, whose processed data it reuses). The
figure is a Bokeh json_item with named models `embedding:source` (x, y and a
colour code per point) and `embedding:mapper`; recolouring sends only the codes
and palette. Hover names are fetched on demand (as on the map), and selections
are exchanged with the map through the layer's native selection.
"""
import threading
import traceback

import numpy as np
import pandas as pd
from bokeh.embed import json_item
from bokeh.models import ColumnDataSource, CustomJSHover, HoverTool, LinearColorMapper
from bokeh.plotting import figure
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel

from ...core import embedding
from ...core.groups import GROUPS_COL
from ...core.encoding import CategoricalColor
from ...settings.config import embedding_kwargs
from ...startup.utils.progress import reporting
from ..figures import THEME, emit
from . import current_session
from .layers import _view

router = APIRouter()
ALIASES = {"bins": "bins", "spots": "spots", "cells": "cells"}


def _state(view):
    state = view.embedding
    clustering = view.clustering
    return dict(method=state["method"], n=0 if state["names"] is None else len(state["names"]),
                version=state["version"], methods=embedding.METHODS, job=view.job,
                default_colour_by=clustering.cluster_col if clustering.is_clustered else "")


def _embedded(view):
    if view.embedding["names"] is None:
        raise HTTPException(404, "Run an embedding first.")
    return view.embedding


def _palette(view, column):
    """Assigned colours of a category column: clusters, groups, or the map's own."""
    if view.clustering.is_clustered and column == view.clustering.cluster_col:
        return dict(view.clustering._color_map)
    if column == GROUPS_COL:
        return dict(view.groups._color_map)
    param = view.rows.get_param("color")
    if isinstance(param, CategoricalColor) and param.colname == column:
        return dict(getattr(param, "category_colors", {}) or {})
    return {}


def _colour_codes(view, column):
    """(palette, codes) of the embedded points; stores the hover values (new version)."""
    state = _embedded(view)
    with view.sample.fetch_lock:
        colours, values = embedding.colours(view.sample, column, state["names"], _palette(view, column))
    codes, palette = pd.factorize(pd.Series(colours, dtype=object))
    state["values"] = values
    state["version"] += 1
    return list(palette), codes.astype(np.uint16)


@router.get("/layers/{layer_id}/embedding")
def state(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        return _state(view)


class EmbeddingRequest(BaseModel):
    method: str = "pca"


@router.post("/layers/{layer_id}/embedding", status_code=202)
def run(layer_id: str, body: EmbeddingRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    if body.method not in embedding.METHODS:
        raise HTTPException(422, f"Unsupported embedding method: {body.method}")
    label = embedding.METHODS[body.method]
    with view.lock:
        if view.job["state"] == "running":
            raise HTTPException(409, "An analysis is already running on this layer.")
        view.job = dict(name="embedding", state="running", message=f"{label} embedding")

    def publish(**changes):
        view.job.update(changes)
        session.publish(dict(type="job", layer=layer_id, **view.job))

    def work():
        try:
            with reporting(lambda message: publish(message=str(message))):
                df = embedding.frame(view.sample, view.clustering, body.method, **embedding_kwargs())
            with view.lock:
                view.embedding.update(method=body.method, names=df["obs_name"].astype(str).to_numpy(),
                                      x=df["x"].to_numpy(np.float32), y=df["y"].to_numpy(np.float32), values=None)
                view.embedding["version"] += 1
            publish(state="done", message=f"Loaded {len(df):,} {ALIASES[view.kind][:-1]} embeddings.")
        except Exception as exc:
            traceback.print_exc()
            publish(state="failed", message=str(exc))

    threading.Thread(target=work, daemon=True, name=f"embedding-{layer_id}").start()
    return view.job


@router.get("/layers/{layer_id}/embedding/figure")
def embedding_figure(layer_id: str, colour_by: str = "", session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        state = _embedded(view)
        palette, codes = _colour_codes(view, colour_by)
        x, y, method, version = state["x"], state["y"], state["method"], state["version"]
    label = embedding.METHODS[method]
    source = ColumnDataSource(dict(x=x, y=y, c=codes), name="embedding:source")
    mapper = LinearColorMapper(palette=palette, low=-.5, high=len(palette) - .5, name="embedding:mapper")
    plot = figure(title=f"{view.title} embedding", height=350, sizing_mode="stretch_width",
                  tools="pan,wheel_zoom,box_select,lasso_select,tap,save,reset", active_scroll="wheel_zoom",
                  output_backend="webgl")
    plot.scatter("x", "y", source=source, size=4, alpha=.8, color={"field": "c", "transform": mapper})
    key = f"embedding:{layer_id}"
    plot.add_tools(HoverTool(name="embedding:hover", tooltips=[("", "$index{custom}")],
                             formatters={"$index": CustomJSHover(
                                 code=f"return window.hdv ? window.hdv.hoverName({key!r}, value) : ''")}))
    plot.xaxis.axis_label, plot.yaxis.axis_label = f"{label} 1", f"{label} 2"
    plot.toolbar.logo = None
    source.selected.js_on_change("indices", emit("embedding-selection", f"{{layer: {layer_id!r}, indices: cb_obj.indices}}"))
    item = json_item(plot, "embedding", theme=THEME)
    item["version"] = version
    return item


@router.get("/layers/{layer_id}/embedding/colours")
def colours(layer_id: str, colour_by: str = Query(""), session=Depends(current_session)):
    """Recolour without resending coordinates: palette and one code per point."""
    view = _view(session, layer_id)
    with view.lock:
        palette, codes = _colour_codes(view, colour_by)
        return dict(version=view.embedding["version"], palette=palette, codes=codes.tolist())


class NamesRequest(BaseModel):
    row_version: int
    rows: list[int]


@router.post("/layers/{layer_id}/embedding/names")
def names(layer_id: str, body: NamesRequest, session=Depends(current_session)):
    """Hover text (native name and colour value) of embedded points."""
    view = _view(session, layer_id)
    with view.lock:
        state = view.embedding
        if state["names"] is None or body.row_version != state["version"]:
            return dict(row_version=state["version"], rows=[])
        names, values = state["names"], state["values"] or [""] * len(state["names"])
        rows = [[i, f"{names[i]}  ·  {values[i]}" if values[i] else str(names[i])]
                for i in body.rows if 0 <= i < len(names)]
        return dict(row_version=state["version"], rows=rows)


class SelectionRequest(BaseModel):
    indices: list[int]


@router.post("/layers/{layer_id}/embedding/selection")
def select(layer_id: str, body: SelectionRequest, session=Depends(current_session)):
    """Points selected on the embedding become the layer's selection (shown on the map)."""
    view = _view(session, layer_id)
    with view.lock:
        names = _embedded(view)["names"]
        view.selected = {str(names[i]) for i in body.indices if 0 <= i < len(names)}
        return dict(row_version=view.rows.row_version, rows=view.selected_positions(), selected=len(view.selected))


@router.get("/layers/{layer_id}/embedding/selection")
def selection(layer_id: str, session=Depends(current_session)):
    """Positions of the layer's selection among the embedded points."""
    view = _view(session, layer_id)
    with view.lock:
        names = _embedded(view)["names"]
        if not view.selected:
            return dict(indices=[])
        return dict(indices=np.flatnonzero(pd.Index(names).isin(list(view.selected))).tolist())
