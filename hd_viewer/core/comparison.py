"""Sample comparison around analysis.comparison.compare_samples, without any UI.

Local metrics are stored on the low-resolution (reference) sample as
`cmp_<metric>` columns so they can be shown on the map like any metadata. The
HTTP API (routers.comparison) uses these functions.
"""
import numpy as np
import pandas as pd

from .. import settings
from ..analysis.comparison import oriented_name

PREFIX = settings.COMPARISON_OBS_PREFIX
GENE_COLUMNS = [
    ("spatial_pearson", "Spatial r"),
    ("capture_ratio", "Capture"),
    ("morans_i_sample1", "Moran low"),
    ("morans_i_sample2", "Moran high"),
    ("zero_fraction_sample1", "Zeros low"),
    ("zero_fraction_sample2", "Zeros high"),
]


def by_area(samples):
    """Keys of `samples` (key -> Sample), largest geometries (the low-res default) first."""
    return sorted(samples, key=lambda key: -samples[key].gdf.geometry.area.mean())


def label_columns(sample):
    """Categorical `adata.obs` columns (cluster labels), without comparison metrics."""
    if sample is None:
        return []
    obs = sample.adata.obs
    return [c for c in obs.columns if not c.startswith(PREFIX) and not pd.api.types.is_numeric_dtype(obs[c])]


def proportion_keys(sample):
    if sample is None:
        return []
    obsm = sample.adata.obsm
    return [k for k in obsm.keys() if isinstance(obsm[k], pd.DataFrame)]


def celltype_options(sample):
    """(key, label): deconvolution proportions (obsm DataFrames) first, then label columns."""
    proportions = proportion_keys(sample)
    labels = [c for c in label_columns(sample) if c not in proportions]
    return [(k, f"{k} (proportions)") for k in proportions] + [(c, c) for c in labels]


def flatten_global(global_metrics, swapped):
    """Rows (metric, value, help) of the global metrics; dict values become `key [sub]` rows."""
    from ..help.comparison import global_help
    rows = []
    for key, value in global_metrics.items():
        items = value.items() if isinstance(value, dict) else [(None, value)]
        for sub, val in items:
            metric = key if sub is None else f"{key} [{sub}]"
            rows.append(dict(metric=metric, value=float(val) if val is not None else np.nan,
                             help=global_help(oriented_name(metric, swapped))))
    return rows


def local_help(metric, swapped):
    from ..help.comparison import LOCAL_HELP
    return LOCAL_HELP.get(oriented_name(metric, swapped), "")


def gene_help():
    from ..help.comparison import GENE_HELP
    return {column: GENE_HELP.get(column, "") for column, _ in GENE_COLUMNS}


def write_obs_columns(sample, local, swapped):
    """Store the local metrics in the low-res `adata.obs` (NaN where unmatched); returns their names."""
    obs = sample.adata.obs
    obs_names = obs.index.astype(str)
    metrics = []
    for col in local.columns:
        # The low-res sample's own clusters are already one of its columns
        if col in ("obs_name", oriented_name("cluster_sample1", swapped)):
            continue
        values = pd.Series(local[col].to_numpy(), index=local["obs_name"].to_numpy())
        if not pd.api.types.is_numeric_dtype(values) or pd.api.types.is_bool_dtype(values):
            values = values.map(lambda v: str(v) if pd.notna(v) else np.nan)
        else:
            values = values.astype(float).replace([np.inf, -np.inf], np.nan)
        obs[PREFIX + col] = values.reindex(obs_names).to_numpy()
        metrics.append(col)
    if hasattr(sample, "clear_multiscale_cache"):
        sample.clear_multiscale_cache()
    return metrics


def clear_obs_columns(sample):
    """Drop the metric columns of a previous comparison."""
    obs = sample.adata.obs
    obs.drop(columns=[c for c in obs.columns if c.startswith(PREFIX)], inplace=True)
    if hasattr(sample, "clear_multiscale_cache"):
        sample.clear_multiscale_cache()
