"""Differential expression between a layer's groups and/or clusters.

Runs as a background job (events `job` with name "de"); the results table is
served in pages (search, sort) and the plots as a Bokeh json_item.
"""
import math
import threading
import traceback

from bokeh.embed import json_item
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field

from ... import settings
from ...analysis.de import DEAnalysis, DEPresenter, build_de_plots
from ...core.de_sets import de_colours, de_inputs, de_options
from ...startup.utils.progress import reporting
from . import current_session
from .layers import _view

router = APIRouter()
ALIASES = {"bins": "Bin", "spots": "Spot", "cells": "Cell"}


def _displayed_factor(view):
    sample = view.sample
    get_factor = getattr(sample, "get_aggregation_factor", None)
    return get_factor(view.rows.effective_resolution_level) if get_factor else 1


def set_options(view):
    """The groups and clusters a comparison can choose (names and colours); call under view.lock."""
    groups, clustering = view.groups, view.clustering
    group_names = de_options("groups", groups, clustering)
    cluster_names = de_options("clusters", groups, clustering)
    clusters = de_colours("clusters", groups, clustering) if cluster_names else {}
    return dict(groups=[dict(name=n, colour=groups._color_map.get(n)) for n in group_names],
                clusters=[dict(name=n, colour=clusters.get(n)) for n in cluster_names])


@router.get("/layers/{layer_id}/de/options")
def options(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        factors = list(getattr(view.sample, "AGGREGATION_FACTORS", ()))
        return dict(**set_options(view), factors=factors, displayed_factor=_displayed_factor(view),
                    ready=view.de["header"] is not None, header=view.de["header"])


class DERequest(BaseModel):
    groups: list[str] = Field(default_factory=list)
    clusters: list[str] = Field(default_factory=list)
    mode: int = 0                 # 0: each vs rest, 1: between the selected sets
    factor: int | None = None     # bins: aggregation level (None: the displayed one)


@router.post("/layers/{layer_id}/de", status_code=202)
def run_de(layer_id: str, body: DERequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    chosen = [("groups", name) for name in body.groups] + [("clusters", name) for name in body.clusters]
    if not chosen:
        raise HTTPException(422, "Select at least one group or cluster to compare.")
    if body.mode == 1 and len(chosen) < 2:
        raise HTTPException(422, '"Between selected" needs at least 2 selected groups or clusters.')
    with view.lock:
        if view.job["state"] == "running":
            raise HTTPException(409, "An analysis is already running on this layer.")
        names, masks, colours, source = de_inputs(view.sample, chosen, view.groups, view.clustering)
        if len(masks) != len(names):
            raise HTTPException(422, "A selected group or cluster no longer exists.")
        factor = body.factor or _displayed_factor(view)
        view.job = dict(name="de", state="running", message="Starting…")
        view.de["generation"] += 1
        view.de["running"] = (source, names)
        generation = view.de["generation"]
        analysis = view.de.setdefault("analysis", DEAnalysis(view.sample))

    def publish(**changes):
        view.job.update(changes)
        session.publish(dict(type="job", layer=layer_id, **view.job))

    def work():
        try:
            with reporting(lambda message: publish(message=str(message))):
                result = analysis.run(masks, body.mode, source=source, factor=factor)
                prepared = DEPresenter(alias=ALIASES[view.kind], sample_type=view.kind, colors=colours).prepare(result)
            with view.lock:
                if view.de["generation"] != generation:
                    publish(state="done", message="DE results discarded because the compared groups or clusters changed.")
                    return
                view.de.update(header=prepared.header, merged=prepared.merged, plot_data=prepared.plot_data,
                               names=list(result.results), n_overlap=result.n_overlap, running=None,
                               compared=(source, names), colours={name: colours.get(name) for name in result.results})
            left_out = f", {result.n_overlap} overlapping observations left out" if result.n_overlap else ""
            unit, units = {"groups": ("group", "groups"), "clusters": ("cluster", "clusters")}.get(
                source, ("set", "groups/clusters"))
            count = len(result.results)
            publish(state="done", message=f"DE done: {count} {unit if count == 1 else units}{left_out}")
        except Exception as exc:
            traceback.print_exc()
            with view.lock:
                if view.de["generation"] == generation:
                    view.de["running"] = None
            publish(state="failed", message=str(exc))

    threading.Thread(target=work, daemon=True, name=f"de-{layer_id}").start()
    return view.job


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


@router.get("/layers/{layer_id}/de/table")
def table(layer_id: str, q: str = "", sort: str = "", desc: bool = False,
          offset: int = Query(0, ge=0), limit: int = Query(settings.DE_TABLE_ROWS, ge=1, le=2000),
          session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        merged = view.de["merged"]
        if merged is None:
            raise HTTPException(404, "Run differential expression first.")
        header, names, colours = view.de["header"], view.de["names"], view.de["colours"]
    rows = merged
    if q:
        rows = rows[rows["Gene"].str.contains(q, case=False, regex=False)]
    if sort in rows.columns:
        rows = rows.sort_values(sort, ascending=not desc, na_position="last")
    page = rows.iloc[offset:offset + limit]
    return dict(header=header, sets=names, colours=[colours.get(name) for name in names], total=len(rows),
                genes=len(merged), pvalue_cutoff=settings.DE_PVALUE_CUTOFF,
                columns=list(merged.columns),
                rows=[[record[0]] + [_number(value) for value in record[1:]]
                      for record in page.itertuples(index=False)])


@router.get("/layers/{layer_id}/de/plots")
def plots(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        plot_data = view.de["plot_data"]
        if plot_data is None:
            raise HTTPException(404, "Run differential expression first.")
    if not plot_data:          # build_de_plots would say so in a Div (a widget the page does not load)
        return dict(message="No positive marker genes to plot.")
    return json_item(build_de_plots(plot_data), "de")
