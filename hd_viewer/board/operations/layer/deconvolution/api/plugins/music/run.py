"""Sparse AnnData-to-R bridge. Only bounded target batches become dense in R."""
import json
import os
from pathlib import Path
import subprocess

import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from vhd_deconv import start


def run(job):
    target, reference = job.read_inputs()
    bridge = job.output_dir / 'r_inputs'
    bridge.mkdir()
    mmwrite(bridge / 'reference.mtx', sparse.coo_matrix(reference.X.T))
    mmwrite(bridge / 'target.mtx', sparse.coo_matrix(target.X.T))
    reference.obs[[job.label_column, job.params['sample_column']]].to_csv(bridge / 'reference_obs.tsv', sep='\t')
    (bridge / 'spec.json').write_text(json.dumps(dict(
        params=job.params, label_column=job.label_column, genes=reference.var_names.tolist(),
        target_genes=target.var_names.tolist(), target_ids=target.obs_names.tolist())))
    job.report(.1, 'Running MuSiC weighted deconvolution in R')
    subprocess.run([os.environ.get('VHD_MUSIC_RSCRIPT', 'Rscript'), str(Path(__file__).with_name('run.R')),
                    str(bridge.resolve())], check=True)
    # R writes generated position identifiers; biological barcodes remain in Python.
    fractions = pd.read_csv(bridge / 'fractions.tsv', sep='\t', index_col=0)
    expected = [f's{i}' for i in range(target.n_obs)]
    if list(fractions.index) != expected:
        raise ValueError('MuSiC output observation order changed.')
    fractions.index = target.obs_names
    genes = json.loads((bridge / 'genes.json').read_text())
    job.write_result(fractions, genes)


if __name__ == '__main__':
    run(start())
