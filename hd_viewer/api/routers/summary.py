"""Dataset summary: per-layer statistics, QC figures and marker genes."""
import math

import numpy as np
from bokeh.embed import json_item
from bokeh.layouts import column, row
from bokeh.plotting import figure
from fastapi import APIRouter, Depends

from ... import settings
from ...core import summary
from ..figures import THEME
from . import current_session
from .layers import _loaded, _view

router = APIRouter()
UNITS = {"bins": "bin", "spots": "spot", "cells": "cell"}


@router.get("/summary")
def dataset_summary(session=Depends(current_session)):
    _loaded(session)
    layers = []
    for view in session.views.values():
        with view.lock, view.sample.fetch_lock:
            summary.ensure_qc(view.sample.adata)          # mean UMIs need total_counts
            layers.append(dict(id=view.id, title=view.title, kind=view.kind, **summary.general_stats(view.sample)))
    return dict(layers=layers)


def _histogram(values, title, x_label, unit, colour):
    hist, edges = np.histogram(values, bins=30)
    plot = figure(title=title, height=190, sizing_mode="stretch_width", x_axis_label=x_label,
                  y_axis_label=f"N {unit}s", tools="reset,pan,wheel_zoom")
    plot.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color=colour, line_color="#1e1e2e", alpha=.8)
    return plot


def _scatter(x, y, title, x_label, y_label, colour):
    plot = figure(title=title, height=190, sizing_mode="stretch_width", x_axis_label=x_label, y_axis_label=y_label,
                  tools="reset,pan,wheel_zoom", output_backend="webgl")
    plot.scatter(x, y, size=3, color=colour, alpha=.6)
    return plot


@router.get("/layers/{layer_id}/qc")
def qc(layer_id: str, session=Depends(current_session)):
    """% MT and genes per observation, and both against UMIs."""
    view = _view(session, layer_id)
    with view.lock, view.sample.fetch_lock:
        obs = summary.ensure_qc(view.sample.adata)
        counts, mt, genes = (obs[name].to_numpy(dtype=float)
                             for name in ("total_counts", "pct_counts_mt", "n_genes_by_counts"))
    unit = UNITS[view.kind]
    plots = [_histogram(mt[np.isfinite(mt)], "% MT", "% MT", unit, "#f14e4e"),
             _histogram(genes[np.isfinite(genes)], f"Genes per {unit}", "N genes", unit, "#4e9af1")]
    shown = np.arange(len(counts))
    if len(shown) > settings.QC_SCATTER_MAX_POINTS:     # a random subset keeps the page light
        shown = np.sort(np.random.default_rng(0).choice(shown, settings.QC_SCATTER_MAX_POINTS, replace=False))
    plots += [_scatter(counts[shown], mt[shown], "UMIs vs % MT", "UMIs", "% MT", "#4ef1a1"),
              _scatter(counts[shown], genes[shown], "UMIs vs Genes", "UMIs", "N genes", "#f1a14e")]
    for plot in plots:
        plot.xaxis.major_label_orientation = .8
        plot.xaxis.ticker.desired_num_ticks = 5
        plot.toolbar.logo = None
    layout = column(row(*plots[:2], sizing_mode="stretch_width"), row(*plots[2:], sizing_mode="stretch_width"),
                    sizing_mode="stretch_width")
    item = json_item(layout, "qc", theme=THEME)
    item["shown"], item["total"] = int(len(shown)), int(len(counts))
    return item


@router.get("/layers/{layer_id}/markers")
def markers(layer_id: str, session=Depends(current_session)):
    """Top marker genes per cluster of the layer's last clustering."""
    view = _view(session, layer_id)
    with view.lock:
        clustering = view.clustering
        found = clustering.get_marker_genes(n_genes=settings.METADATA_MARKER_GENE_COUNT) if clustering.is_clustered else {}
    rows = [dict(cluster=str(cluster), gene=str(gene), score=float(score),
                 pval=None if pval is None or math.isnan(pval) else float(pval), logfc=float(logfc))
            for cluster, genes in found.items() for gene, score, pval, logfc in genes]
    return dict(clustered=clustering.is_clustered, markers=rows)
