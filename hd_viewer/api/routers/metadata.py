"""Descriptive comparison of a layer's numeric metadata (e.g. deconvolution
fractions) between groups and/or clusters, or each against The rest.

Summaries use all native observations; the charts are served as a Bokeh
json_item, kept until the compared groups or clusters change.
"""
import numpy as np
from bokeh.embed import json_item
from bokeh.layouts import column
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field

from ...analysis.metadata_comparison import (build_metadata_plot, compare_metadata,
                                                           compare_metadata_to_rest, metadata_presets,
                                                           numeric_columns, rest_mask)
from ...core.de_sets import de_inputs, de_masks, de_options
from ..figures import THEME
from . import current_session
from .de import set_options
from .layers import _view

router = APIRouter()
NUMERIC = "numeric"
HINT = ("Uses all native observations, independently of the viewport or displayed bin aggregation. "
        "Overlapping observations contribute to each group. Missing and infinite values are omitted per column; "
        "bars show descriptive summaries, without significance tests.")


def _presets(view):
    """Column sets: numeric metadata, then each deconvolution run's readouts."""
    adata = view.sample.adata
    specs = metadata_presets(adata)
    presets = [dict(key=NUMERIC, title="Numeric metadata", unit="value", fractions=False, all_columns=False,
                    columns=[dict(name=name, label=name) for name in numeric_columns(adata.obs)])]
    presets += [dict(key=key, title=spec["title"], unit=spec["unit"], fractions=spec["fractions"], all_columns=True,
                     columns=[dict(name=name, label=label) for name, label in spec["labels"].items()])
                for key, spec in specs.items()]
    # The latest run, preferring its fractions
    default = NUMERIC
    if specs:
        default = next(reversed(specs))
        fraction = default.rsplit(":", 1)[0] + ":fraction"
        default = fraction if fraction in specs else default
    return presets, default


@router.get("/layers/{layer_id}/metadata/options")
def options(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        presets, default = _presets(view)
        return dict(**set_options(view), presets=presets, default_preset=default, hint=HINT,
                    ready=view.metadata["item"] is not None)


class CompareRequest(BaseModel):
    groups: list[str] = Field(default_factory=list)
    clusters: list[str] = Field(default_factory=list)
    include_rest: bool = False
    preset: str = NUMERIC
    columns: list[str] = Field(default_factory=list)
    statistic: str = "mean"


@router.post("/layers/{layer_id}/metadata/compare")
def compare(layer_id: str, body: CompareRequest, session=Depends(current_session)):
    view = _view(session, layer_id)
    chosen = [("groups", name) for name in body.groups] + [("clusters", name) for name in body.clusters]
    if len(chosen) < (1 if body.include_rest else 2) or not body.columns:
        requirement = ("at least one group or cluster with The rest" if body.include_rest
                       else "at least two groups or clusters")
        raise HTTPException(422, f"Select {requirement} and one or more numeric metadata columns.")
    if body.statistic not in ("mean", "median"):
        raise HTTPException(422, "Choose mean or median.")
    sample = view.sample
    with view.lock:
        presets, _ = _presets(view)
        spec = next((preset for preset in presets if preset["key"] == body.preset), None)
        if spec is None:
            raise HTTPException(422, "Unknown metadata set.")
        labels = {column["name"]: column["label"] for column in spec["columns"]}
        if set(body.columns) - set(labels):
            raise HTTPException(422, "Selected metadata columns must exist and contain numeric values.")
        generation = view.metadata["generation"]
        _, masks, colours, _ = de_inputs(sample, chosen, view.groups, view.clustering)
        if len(masks) != len(chosen):
            raise HTTPException(422, "A selected group or cluster no longer exists.")
        with sample.fetch_lock:
            union = np.logical_or.reduce(list(masks.values()))
            defined = None
            if body.include_rest:
                names = de_options("groups", view.groups, view.clustering)
                defined = de_masks(sample, "groups", names, view.groups, view.clustering)
                union = union | rest_mask(defined, sample.adata.n_obs)
            frame = sample.adata.obs.loc[union, body.columns].copy()

    # Summaries outside the lock (millions of native bins on 2 µm data)
    try:
        selected = {name: mask[union] for name, mask in masks.items()}
        if body.include_rest:
            comparisons = compare_metadata_to_rest(frame, selected, {name: mask[union] for name, mask in defined.items()},
                                                   body.columns, body.statistic)
        else:
            comparisons = {"Selected groups": compare_metadata(frame, selected, body.columns, body.statistic)}
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # No Div title (the page loads no BokehJS widgets bundle): the results panel names the layer
    charts = []
    for name, comparison in comparisons.items():
        pair_colours = dict(colours)
        if body.include_rest:
            pair_colours[comparison.groups[-1]] = "#888888"
        chart = build_metadata_plot(comparison, labels if spec["key"] != NUMERIC else None, pair_colours,
                                    spec["unit"], spec["fractions"])
        if body.include_rest:
            chart.title.text = f"{name} vs {comparison.groups[-1]}"
        charts.append(chart)
    item = json_item(column(*charts, sizing_mode="stretch_width"), "metadata", theme=THEME)
    with view.lock:
        if view.metadata["generation"] != generation:      # groups or clusters changed meanwhile
            raise HTTPException(409, "The groups or clusters changed during the comparison; run it again.")
        view.metadata.update(item=item, compared=dict(chosen=chosen, include_rest=body.include_rest))
    return dict(ready=True, plots=len(charts))


@router.get("/layers/{layer_id}/metadata/plots")
def plots(layer_id: str, session=Depends(current_session)):
    view = _view(session, layer_id)
    with view.lock:
        item = view.metadata["item"]
    if item is None:
        raise HTTPException(404, "Run a metadata comparison first.")
    return item
