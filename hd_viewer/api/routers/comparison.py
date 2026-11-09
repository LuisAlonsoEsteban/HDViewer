"""Comparison of two layers (samples) of the same tissue at different resolutions.

The comparison is a background job of the session (events `job` with layer
"workspace" and name "comparison"). Local metrics become `cmp_<metric>` columns
of the low-resolution (reference) layer, shown on the map with its colour-by;
a previous run's columns are replaced. Histogram and scatter are a Bokeh
json_item; selecting scatter points selects those observations on the map.
"""
import math
import threading
import traceback

import numpy as np
import pandas as pd
from bokeh.embed import json_item
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel

from ... import settings
from ...analysis.comparison import compare_samples, oriented_name
from ...core import comparison as core
from ...settings.config import comparison_kwargs
from ...startup.utils.progress import reporting
from ..figures import THEME, emit
from . import current_session
from .layers import _loaded

router = APIRouter()
WORKSPACE = "workspace"          # the `layer` of the comparison's job events
UNITS = {"bins": "bin", "spots": "spot", "cells": "cell"}


def _state(session):
    state = session.comparison
    if not state:
        state.update(job=dict(name="comparison", state="idle", message=""), result=None, low=None, high=None,
                     swapped=False, sample1=None, sample2=None, metrics=[])
    return state


def _comparable(session):
    """Layers with expression data, keyed by id."""
    _loaded(session)
    return {lid: view for lid, view in session.views.items() if view.sample.has_adata()}


@router.get("/comparison")
def options(session=Depends(current_session)):
    views = _comparable(session)
    state = _state(session)
    order = core.by_area({lid: view.sample for lid, view in views.items()})
    layers = []
    for lid in views:
        sample = views[lid].sample
        layers.append(dict(id=lid, title=views[lid].title, clusters=core.label_columns(sample),
                           celltypes=[dict(value=key, label=label) for key, label in core.celltype_options(sample)]))
    return dict(layers=layers, default_sample1=order[0] if order else None,
                default_sample2=order[1] if len(order) > 1 else None,
                predicates=[dict(value="contains", label="inside"), dict(value="intersects", label="overlapping")],
                predicate=settings.COMPARISON_PREDICATE, min_coverage=settings.COMPARISON_MIN_COVERAGE,
                coverage_step=settings.COMPARISON_COVERAGE_STEP, job=state["job"],
                ready=state["result"] is not None)


class CompareRequest(BaseModel):
    sample1: str
    sample2: str
    predicate: str = settings.COMPARISON_PREDICATE
    min_coverage: float = settings.COMPARISON_MIN_COVERAGE
    sample1_cluster_key: str | None = None
    sample2_cluster_key: str | None = None
    sample1_celltype_key: str | None = None
    sample2_celltype_key: str | None = None


def _refresh_metric_layers(session, views):
    """Layers coloured by a comparison column re-read it (or fall back to a gene if it is gone)."""
    for view in views:
        with view.lock:
            name = view.colour_by.get("name") or ""
            if not name.startswith(core.PREFIX):
                continue
            if name in view.sample.get_colnames():
                view.rows.update_colname(name, name)
            else:
                view.colour_by_column(view.sample.adata.var_names[0])
            view.rows.sink.take()
        session.publish(dict(type="layer_stale", layer=view.id))


@router.post("/comparison", status_code=202)
def compare(body: CompareRequest, session=Depends(current_session)):
    views = _comparable(session)
    first, second = views.get(body.sample1), views.get(body.sample2)
    if first is None or second is None or first is second:
        raise HTTPException(422, "Choose two different samples.")
    if body.predicate not in ("contains", "intersects"):
        raise HTTPException(422, "Choose inside or overlapping.")
    if not 0 <= body.min_coverage <= 1:
        raise HTTPException(422, "The minimum coverage must be between 0 and 1.")
    state = _state(session)
    with session.lock:
        if state["job"]["state"] == "running":
            raise HTTPException(409, "A comparison is already running.")
        state["job"] = dict(name="comparison", state="running", message=f"Comparing {first.title} with {second.title}")
    kwargs = dict({k: v for k, v in comparison_kwargs().items() if k not in ("min_coverage", "predicate")},
                  min_coverage=body.min_coverage, predicate=body.predicate,
                  sample1_cluster_key=body.sample1_cluster_key, sample2_cluster_key=body.sample2_cluster_key,
                  sample1_celltype_key=body.sample1_celltype_key, sample2_celltype_key=body.sample2_celltype_key)

    def publish(**changes):
        state["job"].update(changes)
        session.publish(dict(type="job", layer=WORKSPACE, **state["job"]))

    def work():
        try:
            with reporting(lambda message: publish(message=str(message))):
                result = compare_samples(first.sample, second.sample, **kwargs)
            swapped = result.get("reference_sample", 1) == 2
            low, high = (second, first) if swapped else (first, second)
            previous = state["low"]
            touched = {view for view in (previous and views.get(previous), low) if view is not None}
            metrics = []
            for view in touched:                 # results of a previous run are replaced
                with view.lock, view.sample.fetch_lock:
                    core.clear_obs_columns(view.sample)
            if result["local"] is not None:
                with low.lock, low.sample.fetch_lock:
                    metrics = core.write_obs_columns(low.sample, result["local"], swapped)
            state.update(result=result, low=low.id, high=high.id, swapped=swapped, sample1=first.id,
                         sample2=second.id, metrics=metrics)
            _refresh_metric_layers(session, touched)
            low_unit, high_unit = UNITS[low.kind], UNITS[high.kind]
            if result["local"] is None:
                message = (f"No {low_unit} passes the coverage filter: only tissue overlap metrics are available.")
            else:
                matched = int(result["global"][oriented_name("n_sample2_matched", swapped)])
                message = (f"Compared {len(result['local']):,} {low_unit}s with {matched:,} {high_unit}s. Local metrics "
                           f"are available as {core.PREFIX}* columns of the {low_unit}s.")
            publish(state="done", message=message)
        except Exception as exc:
            traceback.print_exc()
            publish(state="failed", message=str(exc))

    threading.Thread(target=work, daemon=True, name=f"comparison-{session.sid[:6]}").start()
    return state["job"]


def _result(session):
    state = _state(session)
    if state["result"] is None:
        raise HTTPException(404, "Compare two samples first.")
    return state


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) or math.isinf(value) else value


def _numeric_metrics(state):
    local = state["result"]["local"]
    return [] if local is None else [m for m in state["metrics"] if pd.api.types.is_numeric_dtype(local[m])]


@router.get("/comparison/result")
def result(session=Depends(current_session)):
    state = _result(session)
    result, swapped = state["result"], state["swapped"]
    numeric = _numeric_metrics(state)
    genes = []
    columns = [dict(key=key, title=title, help=core.gene_help()[key]) for key, title in core.GENE_COLUMNS]
    if result["genes"] is not None:
        table = result["genes"].sort_values(oriented_name("spatial_pearson", swapped), ascending=False)
        values = {key: table[oriented_name(key, swapped)].to_numpy(dtype=float) for key, _ in core.GENE_COLUMNS}
        genes = [[str(gene)] + [_number(values[key][i]) for key, _ in core.GENE_COLUMNS]
                 for i, gene in enumerate(table.index)]
    views = session.views
    x, y = oriented_name("total_counts_sample1", swapped), oriented_name("total_counts_sample2", swapped)
    return dict(sample1=state["sample1"], sample2=state["sample2"], low=state["low"], high=state["high"],
                low_title=views[state["low"]].title, high_title=views[state["high"]].title,
                n_local=0 if result["local"] is None else len(result["local"]),
                global_metrics=[dict(row, value=_number(row["value"]))
                                for row in core.flatten_global(result["global"], swapped)],
                metrics=numeric, metric_help={m: core.local_help(m, swapped) for m in numeric},
                default_metric="pearson_hvg" if "pearson_hvg" in numeric else (numeric[0] if numeric else None),
                default_x=x if x in numeric else None, default_y=y if y in numeric else None,
                gene_columns=columns, genes=genes, prefix=core.PREFIX)


def _local(state, metric):
    local = state["result"]["local"]
    if local is None or metric not in _numeric_metrics(state):
        raise HTTPException(422, f"Unknown local metric: {metric}")
    return pd.to_numeric(local[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _scatter_rows(state, x, y):
    """Local rows with both metrics (the scatter's points, in order)."""
    xs, ys = _local(state, x), _local(state, y)
    ok = (xs.notna() & ys.notna()).to_numpy()
    return xs[ok].to_numpy(), ys[ok].to_numpy(), state["result"]["local"]["obs_name"].to_numpy()[ok]


@router.get("/comparison/plots")
def plots(metric: str, x: str, y: str, session=Depends(current_session)):
    """Histogram of a local metric and a scatter of two of them."""
    state = _result(session)
    unit = UNITS[session.views[state["low"]].kind]
    values = _local(state, metric).dropna()
    histogram = figure(height=200, sizing_mode="stretch_width", tools="pan,wheel_zoom,reset,save",
                       x_axis_label=metric, y_axis_label=f"N {unit}s", title=metric)
    if len(values):
        hist, edges = np.histogram(values, bins=30)
        histogram.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="#4e9af1",
                       line_color="#1e1e2e", alpha=.8)
    xs, ys, names = _scatter_rows(state, x, y)
    source = ColumnDataSource(dict(x=xs, y=ys, obs_name=names.astype(str)), name="comparison:source")
    scatter = figure(height=250, sizing_mode="stretch_width", tools="pan,wheel_zoom,box_select,lasso_select,reset,save",
                     x_axis_label=x, y_axis_label=y, output_backend="webgl")
    scatter.scatter("x", "y", source=source, size=4, alpha=.6, color="#f1a14e")
    scatter.add_tools(HoverTool(tooltips=[("", "@obs_name"), ("x", "@x"), ("y", "@y")]))
    source.selected.js_on_change("indices", emit("comparison-selection", "{indices: cb_obj.indices}"))
    for plot in (histogram, scatter):
        plot.toolbar.logo = None
    return json_item(column(histogram, scatter, sizing_mode="stretch_width"), "comparison", theme=THEME)


class SelectionRequest(BaseModel):
    x: str
    y: str
    indices: list[int]


@router.post("/comparison/selection")
def select(body: SelectionRequest, session=Depends(current_session)):
    """Scatter points selected: those low-res observations become the layer's selection."""
    state = _result(session)
    _, _, names = _scatter_rows(state, body.x, body.y)
    view = session.views[state["low"]]
    with view.lock:
        view.selected = {str(names[i]) for i in body.indices if 0 <= i < len(names)}
        return dict(layer=view.id, row_version=view.rows.row_version, rows=view.selected_positions(),
                    selected=len(view.selected))
