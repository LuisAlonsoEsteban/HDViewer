"""Reference-based NNLS of linear, nonnegative expression (never log counts)."""


import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls

from .progress import report


def _normalise(matrix):
    matrix = sparse.csr_matrix(matrix, dtype=np.float64)
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("Expression must contain finite, nonnegative linear counts.")
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    return sparse.diags(np.divide(1., totals, out=np.zeros_like(totals), where=totals > 0)) @ matrix


def deconvolve(target, reference, label_column, batch_size=256):
    """Return fractions without mutating either input; densify only small batches.

    Mean library-normalised reference profiles are fit on shared genes, then
    NNLS coefficients are normalised to sum to one. These are expression mixture
    estimates, not calibrated cell counts. All-zero observations have NaN fractions.
    """
    if batch_size < 1:
        raise ValueError("Batch size must be positive.")
    if target.n_obs == 0 or reference.n_obs == 0:
        raise ValueError("Target and reference must contain observations.")
    if not target.var_names.is_unique or not reference.var_names.is_unique:
        raise ValueError("Gene names must be unique in both datasets.")
    if label_column not in reference.obs:
        raise ValueError("Choose a cell-type column from the reference.")

    labels = reference.obs[label_column]
    if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
        raise ValueError("Every reference cell must have a cell-type label.")
    labels = labels.astype(str).to_numpy()
    types = sorted(set(labels))
    genes = target.var_names.intersection(reference.var_names, sort=False)

    if len(genes) < max(2, len(types)):
        raise ValueError("Too few shared genes to distinguish the reference cell types.")

    report(f"Building {len(types)} reference profiles from {len(genes)} shared genes")
    
    ref = _normalise(reference[:, genes].X)
    profiles = np.column_stack([
        np.asarray(ref[labels == name].mean(axis=0)).ravel() for name in types
    ])

    if (profiles.sum(axis=0) == 0).any():
        raise ValueError("A reference cell type has no expression on shared genes.")
    if np.linalg.matrix_rank(profiles) < len(types):
        raise ValueError("Reference cell-type profiles are indistinguishable on shared genes.")

    fractions = np.full((target.n_obs, len(types)), np.nan)
    target_columns = target.var_names.get_indexer(genes)
    for start in range(0, target.n_obs, batch_size):
        end = min(start + batch_size, target.n_obs)
        batch = _normalise(target.X[start:end, :][:, target_columns]).toarray()
        for row, expression in enumerate(batch, start):
            if not expression.any():
                continue
            weights, _ = nnls(profiles, expression, maxiter=max(100, 10 * len(types)))
            if weights.sum() > 0:
                fractions[row] = weights / weights.sum()
        report(f"Deconvolved {end:,} / {target.n_obs:,} observations")

    proportions = pd.DataFrame(fractions, index=target.obs_names.copy(), columns=types)
    dominant = pd.Series(np.nan, index=proportions.index, dtype=object)
    assigned = proportions.notna().any(axis=1)
    dominant.loc[assigned] = proportions.loc[assigned].idxmax(axis=1)
    proportions['most_abundant'] = dominant
    return proportions, genes.tolist()
