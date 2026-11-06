"""Two-stage cell2location workflow; dependencies live in this plugin environment."""
import os
import torch
import numpy as np
import pandas as pd
import scvi
from cell2location.models import RegressionModel, Cell2location
from vhd_deconv import start


def run(job):
    target, reference = job.read_inputs()
    p = job.params
    scvi.settings.seed = p['seed']
    scvi.settings.num_threads = int(os.environ.get('VHD_DECONV_THREADS', '8'))
    accelerator = 'gpu' if os.environ.get('VHD_DECONV_GPU') == '1' and torch.cuda.is_available() else 'cpu'
    cached = job.reference_cache / 'signatures.json'
    from vhd_deconv import read_table, write_table
    if cached.exists():
        job.report(.05, 'Loaded cached reference signatures')
        signatures = read_table(cached)
    else:
        # Drop genes absent from the reference before fitting the count model.
        reference = reference[:, np.asarray(reference.X.sum(axis=0)).ravel() > 0].copy()
        RegressionModel.setup_anndata(reference, batch_key=p['batch_column'], labels_key=job.label_column)
        model = RegressionModel(reference)
        job.report(.05, 'Training reference regression')
        model.train(max_epochs=p['reference_epochs'], batch_size=p['batch_size'], accelerator=accelerator, device="auto")
        reference = model.export_posterior(reference, sample_kwargs=dict(
            num_samples=p['posterior_samples'], batch_size=p['batch_size'], accelerator=accelerator, device="auto", summary_frequency=1))
        names = reference.uns['mod']['factor_names']
        columns = [f'means_per_cluster_mu_fg_{name}' for name in names]
        if 'means_per_cluster_mu_fg' in reference.varm:
            signatures = reference.varm['means_per_cluster_mu_fg'][columns].copy()
        else:
            signatures = reference.var[columns].copy()
        signatures.columns = list(names)
        write_table(signatures, job.reference_output / 'signatures.json')
    genes = target.var_names.intersection(signatures.index, sort=False)
    if len(genes) < max(2, signatures.shape[1]):
        raise ValueError('Too few shared expressed genes for cell2location.')
    target = target[:, genes].copy()
    original_index = target.obs_names.copy()
    # Count model registration and inference only on observations with shared counts.
    active = np.asarray(target.X.sum(axis=1)).ravel() > 0
    target = target[active].copy()
    if not target.n_obs:
        raise ValueError('No observations have counts on fitted reference genes.')
    Cell2location.setup_anndata(target)
    model = Cell2location(target, cell_state_df=signatures.loc[genes],
                          N_cells_per_location=p['N_cells_per_location'], detection_alpha=p['detection_alpha'])
    job.report(.4, 'Training spatial abundance model')
    model.train(max_epochs=p['spatial_epochs'], batch_size=p['batch_size'], train_size=1,
                accelerator=accelerator, device="auto")
    target = model.export_posterior(target, sample_kwargs=dict(
        num_samples=p['posterior_samples'], batch_size=p['batch_size'], accelerator=accelerator, device="auto", summary_frequency=1))
    abundance = pd.DataFrame(np.asarray(target.obsm['means_cell_abundance_w_sf']),
                             index=target.obs_names, columns=signatures.columns)
    fractions = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0)
    job.write_result(fractions.reindex(original_index), genes, abundance=abundance.reindex(original_index),
                     abundance_units='expected cells per location (posterior mean)')


if __name__ == '__main__':
    run(start())
