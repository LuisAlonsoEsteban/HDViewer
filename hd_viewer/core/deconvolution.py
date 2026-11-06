"""Deconvolution of a layer, without any UI.

Runs go through the shared DeconvolutionService (one queue per process); their
results are published into the sample as numbered runs (`deconv_<run>:<type>`
columns, `adata.uns['deconvolution_runs']`). The Bokeh UI
(board.layout.layer.deconvolution) and the HTTP API both use these functions.
"""
import json
from functools import lru_cache
from urllib.parse import quote

import anndata as ad
import numpy as np
import pandas as pd

from .. import settings

PREFIX = settings.DECONV_RESULT_PREFIX


@lru_cache(maxsize=1)
def default_service():
    """One queue and runner per viewer process, shared by all layers and sessions."""
    from ..board.operations.layer.deconvolution import DeconvolutionService
    from ..settings.config import deconvolution_service_kwargs
    return DeconvolutionService(**deconvolution_service_kwargs())


def runs(sample):
    """Saved runs of a sample: run id -> metadata record."""
    return sample.adata.uns.get('deconvolution_runs', {})


def cache_namespace(sample, kind, key, source_bin_size_um=None):
    """Native input resolution and layer, independent of display aggregation."""
    size = source_bin_size_um or getattr(sample, 'bin_size', None)
    resolution = 'native'
    if size is not None:
        size = float(size)
        if np.isfinite(size) and size > 0:
            label = f'{int(size):03d}' if size.is_integer() else f'{size:g}'.replace('.', 'p')
            resolution = f'square_{label}um'
    # Encode identifiers as individual directory names, including dots/slashes.
    parts = [resolution, quote(str(kind), safe='').replace('.', '%2E')]
    if key != kind:
        parts.append(quote(str(key), safe='').replace('.', '%2E'))
    return '/'.join(parts)


def snapshot_target(sample, method):
    """Counts (and spatial coordinates for RCTD) copied under the sample lock, so
    later changes cannot alter what a background run consumes."""
    with sample.fetch_lock:
        target = ad.AnnData(sample.adata.X.copy(), obs=pd.DataFrame(index=sample.adata.obs_names.copy()),
                            var=pd.DataFrame(index=sample.adata.var_names.copy()))
        if method == 'spacexr' and 'spatial' in sample.adata.obsm:
            target.obsm['spatial'] = sample.adata.obsm['spatial'].copy()
    return target


def publish_result(sample, result, label_column, reference_path='', reference_name=''):
    """Validate `result` and store it in the sample as a new run; returns the run id."""
    from ..board.operations.layer.deconvolution import DeconvolutionResult, validate_result
    if not isinstance(result, DeconvolutionResult):
        fractions, genes = result
        result = DeconvolutionResult(fractions.drop(columns='most_abundant', errors='ignore'), genes,
                                     metadata=dict(method='nnls', version='0.1.0', params={}))

    validate_result(result, sample.adata.obs_names)
    fractions = result.fractions.copy()
    dominant = pd.Series(np.nan, index=fractions.index, dtype=object)
    assigned = fractions.notna().any(axis=1)
    dominant.loc[assigned] = fractions.loc[assigned].idxmax(axis=1)
    fractions['most_abundant'] = dominant

    # Validate everything before publishing under the sample lock.
    with sample.fetch_lock:
        history = runs(sample)
        run = 1
        while str(run) in history or any(str(c).startswith(f'{PREFIX}{run}:') for c in sample.get_colnames()):
            run += 1
        prefix = f'{PREFIX}{run}:'
        fraction_columns = [prefix + name for name in result.fractions.columns]
        abundance_columns = []
        for name, original in zip(fraction_columns, result.fractions.columns):
            sample.adata.obs[name] = result.fractions[original].astype(float)
        sample.adata.obs[prefix + 'most_abundant'] = dominant.astype('category')
        if result.abundance is not None:
            for i, original in enumerate(result.abundance.columns):
                # Separate namespace avoids collisions with biological cell-type names.
                name = f'{PREFIX}{run}_abundance:{i}'
                abundance_columns.append(name)
                sample.adata.obs[name] = result.abundance[original].astype(float)
            sample.adata.obsm[f'deconvolution_{run}_abundance'] = result.abundance.copy()
            sample.adata.obsm['deconvolution_abundance'] = result.abundance.copy()
        else:
            sample.adata.obsm.pop('deconvolution_abundance', None)
        metadata = dict(
            method=result.metadata.get('method', 'nnls'), version=result.metadata.get('version', ''),
            params_json=json.dumps(result.metadata.get('params', {}), sort_keys=True),
            reference_id=result.metadata.get('reference_id', ''), reference_path=reference_path or '',
            reference_name=reference_name or '', label_column=label_column,
            cell_types=list(result.fractions.columns), genes=result.genes,
            fraction_columns=fraction_columns, abundance_columns=abundance_columns,
            abundance_units=result.metadata.get('abundance_units', ''),
            fraction_units=result.metadata.get('fraction_units', 'relative fraction'),
            log_path=result.metadata.get('log_path', ''), run_id=result.metadata.get('run_id', ''),
            cache_namespace=result.metadata.get('cache_namespace', ''),
        )
        sample.adata.obsm['deconvolution'] = fractions
        sample.adata.obsm[f'deconvolution_{run}'] = fractions.copy()
        sample.adata.uns['deconvolution'] = metadata
        sample.adata.uns.setdefault('deconvolution_runs', {})[str(run)] = metadata.copy()
        means = list(sample.adata.uns.get('mean_obs_columns', []))
        sample.adata.uns['mean_obs_columns'] = list(dict.fromkeys(means + fraction_columns))
        if hasattr(sample, 'clear_multiscale_cache'):
            sample.clear_multiscale_cache()
    return str(run)
