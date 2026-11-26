"""Export of a layer's data as an AnnData file (.h5ad), without any UI.

The file holds the layer's expression matrix and every per-observation result
stored in it (clusters, groups, deconvolution proportions, embeddings).
"""
import uuid
from pathlib import Path

import anndata as ad
import pandas as pd


def _writable_obs(obs):
    """`obs` with columns h5ad can store: mixed or object columns become strings
    (missing values stay missing as categories cannot hold None)."""
    obs = obs.copy()
    for column in obs.columns:
        series = obs[column]
        if series.dtype == object:
            obs[column] = pd.Categorical(series.where(series.isna(), series.astype(str)))
    obs.columns = obs.columns.astype(str)
    return obs


_SERVER_KEYS = {"log_path", "cache_namespace"}


def _without_server_paths(value):
    """Metadata without entries that only name files on the server."""
    if isinstance(value, dict):
        return {key: _without_server_paths(item) for key, item in value.items() if key not in _SERVER_KEYS}
    return value


def export_layer(sample, directory):
    """Write `sample.adata` to a new, randomly named file in `directory`; returns its path.

    A lightweight AnnData shares the matrix (no copy). Unstructured metadata
    (`uns`) is kept when h5ad can store it.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{uuid.uuid4().hex}.h5ad"
    with sample.fetch_lock:
        adata = sample.adata
        exported = ad.AnnData(X=adata.X, obs=_writable_obs(adata.obs), var=adata.var.copy(),
                              obsm=dict(adata.obsm), layers=dict(adata.layers))
        exported.uns = _without_server_paths(dict(adata.uns))
        try:
            try:
                exported.write_h5ad(path)
            except Exception:
                exported.uns = {}
                exported.write_h5ad(path)
        except Exception:
            path.unlink(missing_ok=True)
            raise
    return path
