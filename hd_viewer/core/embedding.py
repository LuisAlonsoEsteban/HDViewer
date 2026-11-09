"""2-D embeddings (PCA, UMAP, t-SNE) of a layer, without any UI.

Embeddings reuse the clustering state's processed data and its PCA/neighbour
caches; results are stored in the sample (`adata.obsm`, aligned with the native
observations). The HTTP API (routers.embedding) uses these functions.
"""
from functools import lru_cache

import numpy as np
import pandas as pd

from .. import settings
from ..startup.utils.progress import report

METHODS = {"pca": "PCA", "umap": "UMAP", "tsne": "t-SNE"}
KEYS = {"pca": "X_pca_2d", "umap": "X_umap", "tsne": "X_tsne"}
DEFAULT_COLOUR = "#6baed6"
MISSING_COLOUR = "#555555"


@lru_cache(maxsize=1)
def viridis():
    """The 256-colour viridis ramp (as Bokeh's Viridis256; core does not import Bokeh)."""
    from matplotlib import colormaps
    from matplotlib.colors import to_hex
    return [to_hex(colour) for colour in colormaps["viridis"](np.linspace(0, 1, 256))]


def _store(sample, clustering, method, values):
    processed = clustering.adata_processed
    if processed is None:
        return
    full = np.full((sample.adata.n_obs, 2), np.nan, dtype=float)
    rows = pd.Index(sample.adata.obs_names.astype(str)).get_indexer(pd.Index(processed.obs_names.astype(str)))
    valid = rows >= 0
    full[rows[valid], :] = values[valid, :]
    sample.adata.obsm[KEYS[method]] = full


def compute(sample, clustering, method, n_pcs=30, n_neighbors=15):
    """The (n_obs, 2) embedding of `sample`, computed once per parameters."""
    import scanpy as sc

    method = method.lower()
    if method not in KEYS:
        raise ValueError(f"Unsupported embedding method: {method}")
    key = KEYS[method]
    step = "pca_embedding" if method == "pca" else method
    params = (n_pcs, n_neighbors) if method == "umap" else n_pcs
    if clustering.adata_processed is not None and clustering._is_cached(step, params) and key in sample.adata.obsm:
        report(f"Using cached {method.upper()} embedding")
        return sample.adata.obsm[key]

    if method == "pca":
        adata = clustering._ensure_pca(n_pcs)
        values = adata.obsm["X_pca"][:, :2]
    elif method == "umap":
        adata = clustering._ensure_neighbors(n_pcs, n_neighbors)
        report("Running UMAP...")
        sc.tl.umap(adata)
        values = adata.obsm["X_umap"][:, :2]
    else:
        adata = clustering._ensure_pca(n_pcs)
        report("Running t-SNE...")
        sc.tl.tsne(adata, n_pcs=n_pcs)
        values = adata.obsm["X_tsne"][:, :2]

    _store(sample, clustering, method, np.asarray(values, dtype=float))
    clustering._mark_computed(step, params)
    return sample.adata.obsm[key]


def frame(sample, clustering, method, n_pcs=30, n_neighbors=15):
    """Embedded observations: columns obs_name, x, y (observations without coordinates left out)."""
    values = np.asarray(compute(sample, clustering, method, n_pcs=n_pcs, n_neighbors=n_neighbors), dtype=float)
    df = pd.DataFrame(values, columns=["x", "y"], index=sample.adata.obs_names.astype(str))
    df.index.name = "obs_name"
    return df.reset_index().dropna(subset=["x", "y"])


def colours(sample, colname, obs_names, palette=None):
    """(colours, hover values) of `obs_names` by `colname`, over every embedded
    observation rather than only the map's loaded viewport. Genes follow the
    sample's expression scale; categories use `palette` (label -> colour) first."""
    n = len(obs_names)
    if not colname or colname not in sample.get_colnames():
        return [DEFAULT_COLOUR] * n, [""] * n
    adata = sample.adata
    rows = adata.obs_names.astype(str).get_indexer(pd.Index(obs_names, dtype=str))
    found = rows >= 0
    safe = np.maximum(rows, 0)

    if sample.is_continous(colname):
        if colname in adata.var_names:
            column = adata.X[:, [adata.var_names.get_loc(colname)]]
            counts = column.toarray().ravel() if hasattr(column, "toarray") else np.asarray(column).ravel()
            values = sample.scale_expression(counts[safe], sample.library_sizes()[safe])
        else:
            values = pd.to_numeric(adata.obs[colname], errors="coerce").to_numpy(dtype=float)[safe]
        values = np.where(found, np.asarray(values, dtype=float), np.nan)
        finite = np.isfinite(values)
        low, high = (np.min(values[finite]), np.max(values[finite])) if finite.any() else (0.0, 1.0)
        span = high - low if high > low else 1.0
        positions = np.clip(((np.nan_to_num(values, nan=low) - low) / span * 255).astype(int), 0, 255)
        ramp = viridis()
        result = [ramp[p] if ok else MISSING_COLOUR for p, ok in zip(positions, finite)]
        return result, [f"{v:.3g}" if ok else "NA" for v, ok in zip(values, finite)]

    labels = adata.obs[colname].astype(object).to_numpy()[safe]
    labels = ["NA" if not ok or pd.isna(label) or str(label) in ("None", "nan", "") else str(label)
              for label, ok in zip(labels, found)]
    palette = dict(palette or {})
    for index, label in enumerate(sorted(set(labels) - set(palette) - {"NA"})):
        palette[label] = settings.CATEGORICAL_PALETTE[index % len(settings.CATEGORICAL_PALETTE)]
    return [palette.get(label, MISSING_COLOUR) for label in labels], labels
