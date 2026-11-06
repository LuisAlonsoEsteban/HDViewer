"""Project canonical native-observation selections onto displayed bin aggregates."""

import numpy as np
import pandas as pd

def _summary_factors(names, sample):
    allowed = getattr(sample, "AGGREGATION_FACTORS", ())
    prefix = getattr(sample, "aggregate_obs_prefix", None)
    if not allowed or prefix is None:          # spots and cells are never aggregated
        return set()
    return {int(name.rsplit(":", 1)[1]) for name in names
            if name.startswith(prefix) and name.rsplit(":", 1)[1].isdigit()
            and int(name.rsplit(":", 1)[1]) in allowed}


def aggregation_names(sample, factor):
    """Native obs -> display cell, without aggregating the expression matrix."""
    return sample.aggregation_names(factor)


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
