"""Host-side AnnData bridge for the unmodified official CIBERSORTx image."""
from dataclasses import replace
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy import sparse

from .api.sdk.vhd_deconv import Job


# The runner appends these from the server environment, keeping them out of argv
CREDENTIALS = {'--username': 'VHD_CIBERSORTX_USERNAME', '--token': 'VHD_CIBERSORTX_TOKEN'}


def _write_expression(adata, columns, path):
    """Stream genes in bounded chunks; never densify the single-cell reference."""
    with path.open('w') as handle:
        handle.write('Gene\t' + '\t'.join(columns) + '\n')
        for start in range(0, adata.n_vars, 32):
            matrix = adata.X[:, start:start + 32]
            matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
            pd.DataFrame(matrix.T, index=adata.var_names[start:start + 32]).to_csv(
                handle, sep='\t', header=False)


def _read(path):
    return pd.read_csv(path, sep='\t', index_col=0, dtype={0: str}, keep_default_na=False)


def run(plugin, runner, input_dir, output_dir, progress, log_path):
    job = Job(input_dir, output_dir)
    target, reference = job.read_inputs()
    types = sorted(set(reference.obs[job.label_column].astype(str)))
    encoded = {name: f'type{i}' for i, name in enumerate(types)}
    ids = [f's{i}' for i in range(target.n_obs)]
    # CIBERSORTx accepts duplicate column headers denoting each cell's phenotype.
    _write_expression(reference, [encoded[x] for x in reference.obs[job.label_column]], input_dir / 'reference.tsv')
    _write_expression(target, ids, input_dir / 'mixture.tsv')
    p = job.params
    args = ['--mixture', 'mixture.tsv', '--perm', str(p['perm']), '--QN', str(p['quantile_normalization']).upper()]
    cached = job.reference_cache / 'signature.tsv'
    if cached.exists():
        args += ['--sigmatrix', 'reference_artifacts/signature.tsv']
    else:
        # The official single-cell mode builds a signature before estimating fractions.
        args += ['--single_cell', 'TRUE', '--refsample', 'reference.tsv',
                 '--G.min', str(p['g_min']), '--G.max', str(p['g_max']), '--q.value', str(p['q_value']),
                 '--replicates', str(p['replicates']), '--sampling', str(p['sampling']), '--fraction', str(p['fraction'])]
    if p['batch_correction'] != 'none':
        args += [f"--rmbatch{p['batch_correction']}mode", 'TRUE']
        if cached.exists():
            args += ['--refsample', 'reference.tsv']
            sources = list(job.reference_cache.glob('*sourceGEP*'))
            if sources:
                args += ['--sourceGEPs', 'reference_artifacts/' + sources[0].name]
    progress(.1, 'Starting official CIBERSORTx single-cell fractions workflow')
    runner.run(replace(plugin, command=args, secret_args=CREDENTIALS), input_dir, output_dir, progress, log_path)
    candidates, signatures = [], []
    for path in sorted(output_dir.glob('*.txt')) + sorted(output_dir.glob('*.tsv')):
        try:
            frame = _read(path)
        except (ValueError, pd.errors.ParserError):
            continue
        if set(encoded.values()).issubset(frame.columns):
            if set(frame.index) == set(ids):
                candidates.append((path, frame))
            elif set(frame.index).issubset(set(reference.var_names)):
                signatures.append((path, frame))
    if len(candidates) != 1:
        raise ValueError('Expected one CIBERSORTx fraction result with all target observations and cell types.')
    _, frame = candidates[0]
    fractions = frame.loc[ids, list(encoded.values())].astype(float)
    fractions.index, fractions.columns = target.obs_names, types
    if cached.exists():
        signature = _read(cached)
    else:
        if len(signatures) != 1:
            raise ValueError('Expected one generated CIBERSORTx signature matrix.')
        signature_path, signature = signatures[0]
        shutil.copy(signature_path, job.reference_output / 'signature.tsv')
        for path in output_dir.glob('*sourceGEP*'):
            if path.is_file():
                shutil.copy(path, job.reference_output / path.name)
    genes = target.var_names.intersection(signature.index, sort=False).tolist()
    job.write_result(fractions, genes)
