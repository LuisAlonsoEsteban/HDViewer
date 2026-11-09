"""Shared test fixtures: small samples and the API's per-layer view.

`LayerView` (hd_viewer.api.views) is what the interface drives: it holds the
layer's rows (core.layer_state.LayerRows with a recording sink), its clustering
(core.clustering.ClusteringState), groups (core.groups.GroupsState) and ROIs
(core.roi.ROIGeometry). Tests of map, analysis and selection logic build one with
`make_view` and drive it like the HTTP API does.
"""
import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.api.views import LayerView
from hd_viewer.startup.sample import Bins, Cells, Spots

COLOUR = "#6baed6"


def grid_frames(n=4, size=1.0, prefix="b", not_sequenced=0):
    """(gdf, names) of an n x n grid of square geometries, plus `not_sequenced`
    extra small geometries with no expression."""
    names = [f"{prefix}{x}_{y}" for y in range(n) for x in range(n)]
    extra = [f"not_sequenced{i}" for i in range(not_sequenced)]
    polygons = [box(x * size, y * size, (x + 1) * size, (y + 1) * size) for y in range(n) for x in range(n)]
    polygons += [box(1.2 + i, 1.2, 1.4 + i, 1.4) for i in range(not_sequenced)]
    gdf = gpd.GeoDataFrame({"index": names + extra}, geometry=polygons,
                           index=pd.Index(names + extra, name="barcode"))
    return gdf, names


def expression(names, n_genes=2, X=None):
    X = np.arange(len(names) * n_genes, dtype=float).reshape(len(names), n_genes) + 1 if X is None else X
    return ad.AnnData(X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=[f"g{i}" for i in range(X.shape[1])]))


def make_sample(kind="bins", n=4, n_genes=2, not_sequenced=0, size=1.0, X=None, **options):
    """A small Bins/Spots/Cells sample on an n x n grid (unit squares by default)."""
    gdf, names = grid_frames(n, size, prefix={"bins": "b", "spots": "s", "cells": "c"}[kind],
                             not_sequenced=not_sequenced)
    adata = expression(names, n_genes, X)
    if kind == "bins":
        options.setdefault("pixels_per_bin", 1)
        options.setdefault("bin_size", 2)
        return Bins(gdf, adata, **options)
    if kind == "spots":
        return Spots(gdf, adata, **options)
    return Cells(gdf, adata, **options)


def make_view(sample=None, kind="bins", layer_id=None, title=None, colour=COLOUR, **sample_options):
    """The API view of `sample` (default: `make_sample(kind, **sample_options)`)."""
    if sample is None:
        sample = make_sample(kind, **sample_options)
    return LayerView(layer_id or kind, kind, sample, title or kind.title(), colour)


def show(view, x0, x1, y0, y1, width=800, height=800):
    """Load the rows of a view box, as POST /layers/{id}/viewport does; returns the ops."""
    return view.update_viewport(x0, x1, y0, y1, width, height)


def loaded(view):
    """Native or display names of the rows currently loaded in the view."""
    return list(view.rows.loaded_obs_names)


def column(view, name):
    """A column of the view's data source as a NumPy array (or list for objects)."""
    values = view.rows.sink.data[name]
    return np.asarray(values) if isinstance(values, np.ndarray) else list(values)
