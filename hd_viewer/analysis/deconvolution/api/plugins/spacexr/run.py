"""Sparse AnnData bridge to dmcable/spacexr's Reference/SpatialRNA RCTD API."""
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from vhd_deconv import start


def read_fractions(result, target_ids, type_ids, mode):
    """Map mode-specific RCTD predictions onto the complete target/type axes."""
    barcodes = result['barcodes']
    if len(set(barcodes)) != len(barcodes) or set(barcodes) - set(target_ids):
        raise ValueError('RCTD returned duplicate or unknown observation identifiers.')
    fractions = pd.DataFrame(np.nan, index=target_ids, columns=type_ids)

    def assign(barcode, types, weights):
        if len(types) != len(weights) or len(set(types)) != len(types) or set(types) - set(type_ids):
            raise ValueError('RCTD returned invalid cell-type identifiers or weights.')
        values = np.asarray(weights, dtype=float)
        if not len(values) or not np.isfinite(values).all() or (values < 0).any() or values.sum() <= 0:
            raise ValueError('RCTD returned invalid weights for an assigned observation.')
        fractions.loc[barcode] = 0.0
        fractions.loc[barcode, types] = values / values.sum()

    predictions = result['predictions']
    if len(predictions) != len(barcodes):
        raise ValueError('RCTD predictions do not match its observation identifiers.')
    for barcode, prediction in zip(barcodes, predictions):
        if mode == 'doublet':
            classification = prediction['spot_class']
            if classification in ('reject', 'doublet_uncertain'):
                # An uncertain second identity does not define confident fractions.
                continue
            if classification == 'singlet':
                assign(barcode, [prediction['first_type']], [1.0])
            elif classification == 'doublet_certain':
                assign(barcode, [prediction['first_type'], prediction['second_type']], prediction['weights'])
            else:
                raise ValueError(f'Unknown RCTD spot classification: {classification}')
        else:
            assign(barcode, prediction['types'], prediction['weights'])
    return fractions


def run(job):
    target, reference = job.read_inputs()
    for name, data in [('Target', target), ('Reference', reference)]:
        values = data.X.data if sparse.issparse(data.X) else np.asarray(data.X)
        if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.round(values), rtol=0, atol=1e-6):
            raise ValueError(f'{name}: RCTD requires finite, nonnegative raw integer counts in .X.')
    labels = reference.obs[job.label_column].astype(str)
    type_names = sorted(set(labels))

    if len(type_names) < 2:
        raise ValueError('RCTD requires at least two reference cell types.')
    if job.params['UMI_min'] > job.params['UMI_max']:
        raise ValueError('UMI_min must not exceed UMI_max.')
    if job.params['n_max_cells'] < job.params['CELL_MIN_INSTANCE']:
        raise ValueError('n_max_cells must be at least CELL_MIN_INSTANCE.')

    # Safe names avoid R's restrictions on type names and preserve labels such as
    # "NA", numeric strings, slashes, or embedded whitespace exactly on return.
    type_ids = [f'type{i}' for i in range(len(type_names))]
    type_map = dict(zip(type_names, type_ids))
    target_ids = [f's{i}' for i in range(target.n_obs)]
    coordinates = None
    if 'spatial' in target.obsm:
        coordinates = np.asarray(target.obsm['spatial'], dtype=float)
        if coordinates.shape != (target.n_obs, 2) or not np.isfinite(coordinates).all():
            raise ValueError('RCTD spatial coordinates must be a finite observations-by-2 matrix.')
        coordinates = coordinates.tolist()

    bridge = job.output_dir / 'r_inputs'
    bridge.mkdir()
    job.report(.02, 'Exporting sparse RCTD reference and target counts')
    mmwrite(bridge / 'reference.mtx', sparse.coo_matrix(reference.X.T))
    mmwrite(bridge / 'target.mtx', sparse.coo_matrix(target.X.T))
    budget = max(1, int(os.environ.get('VHD_DECONV_THREADS', job.params['max_cores'])))
    params = dict(job.params, max_cores=min(job.params['max_cores'], budget))
    (bridge / 'spec.json').write_text(json.dumps(dict(
        params=params, genes=reference.var_names.tolist(), target_genes=target.var_names.tolist(),
        target_ids=target_ids, reference_ids=[f'r{i}' for i in range(reference.n_obs)],
        reference_labels=labels.map(type_map).tolist(), type_ids=type_ids, coordinates=coordinates,
    ), allow_nan=False))
    (job.output_dir / 'cell_type_labels.json').write_text(json.dumps(dict(zip(type_ids, type_names))))
    (job.output_dir / 'observation_labels.json').write_text(json.dumps(dict(zip(target_ids, target.obs_names))))
    if coordinates is None:
        job.report(.05, 'No spatial coordinates supplied; using RCTD placeholder coordinates (not used by the fit)')
    job.report(.1, 'Running whole-layer spacexr RCTD in R')
    executable = os.environ.get('VHD_SPACEXR_RSCRIPT', 'Rscript')
    
    try:
        subprocess.run([executable, '--vanilla', str(Path(__file__).with_name('run.R')),
                        str(bridge.resolve())], check=True)
    except FileNotFoundError as exc:
        raise RuntimeError('RCTD needs Rscript with dmcable/spacexr, Matrix, and jsonlite. '
                           'Install them or set VHD_SPACEXR_RSCRIPT to that Rscript executable.') from exc

    result = json.loads((bridge / 'rctd.json').read_text())
    fractions = read_fractions(result, target_ids, type_ids, params['doublet_mode'])
    fractions.index = target.obs_names
    fractions.columns = type_names
    
    genes = result['genes']
    if not isinstance(genes, list) or len(set(genes)) != len(genes) or set(genes) - set(target.var_names.intersection(reference.var_names)):
        raise ValueError('RCTD must return unique shared gene identifiers.')

    job.write_result(fractions, genes)


if __name__ == '__main__':
    run(start())
