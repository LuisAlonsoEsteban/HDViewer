"""Project canonical native-observation selections onto displayed bin aggregates."""

import numpy as np
import pandas as pd


def _summary_factors(names, sample):
    allowed = getattr(sample, "AGGREGATION_FACTORS", ())
    return {int(name.rsplit(":", 1)[1]) for name in names
            if name.startswith("summary:") and name.rsplit(":", 1)[1].isdigit()
            and int(name.rsplit(":", 1)[1]) in allowed}


def aggregation_names(sample, factor):
    """Native obs -> display cell, without aggregating the expression matrix."""
    existing = getattr(sample, "_aggregated_cell_names", {}).get(factor)
    if existing is not None:
        return existing
    cache = getattr(sample, "_selection_cell_names", None)
    identity = (id(sample.gdf), id(sample.adata))
    if cache is None or cache[0] != identity:
        sample._selection_cell_names = cache = (identity, {})
    if factor not in cache[1]:
        coordinates = sample.get_coordinates(sample.gdf)
        valid = np.isfinite(coordinates.x) & np.isfinite(coordinates.y)
        labels = pd.Series(None, index=sample.gdf.index, dtype=object)
        x = np.floor(coordinates.loc[valid, "x"] / (sample.native_width * factor)).astype(np.int64)
        y = np.floor(coordinates.loc[valid, "y"] / (sample.native_height * factor)).astype(np.int64)
        labels.loc[valid] = "summary:" + x.astype(str) + ":" + y.astype(str) + f":{factor}"
        labels.index = sample.gdf[sample.on].astype(str)
        cache[1][factor] = labels[~labels.index.duplicated()].reindex(sample.adata.obs_names)
    return cache[1][factor]


def native_selection(sample, names):
    names = {str(name) for name in names}
    if not names:
        return set()
    native = set(sample.adata.obs_names.intersection(list(names)))
    for factor in _summary_factors(names, sample):
        labels = aggregation_names(sample, factor)
        native.update(sample.adata.obs_names[labels.isin(names)])
    return native


def selection_indices(sample, selected, displayed_names):
    displayed = np.asarray(displayed_names, dtype=str)
    if not selected or not len(displayed):
        return []
    projected = set(selected)
    for factor in _summary_factors(displayed, sample):
        labels = aggregation_names(sample, factor)
        projected.update(labels.loc[labels.index.isin(selected)].dropna())
    return np.flatnonzero(np.isin(displayed, list(projected))).tolist()
