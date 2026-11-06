"""RCTD bridge, mode conversion, and dataset-cache contracts without an R install."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from hd_viewer.analysis.deconvolution import DeconvolutionService, discover, validate_result
from hd_viewer.analysis.deconvolution.api.sdk import vhd_deconv as sdk


PLUGIN = discover()['spacexr']
module_spec = importlib.util.spec_from_file_location('spacexr_adapter', PLUGIN.path / 'run.py')
adapter = importlib.util.module_from_spec(module_spec)
with patch.dict('sys.modules', vhd_deconv=sdk):
    module_spec.loader.exec_module(adapter)


class SpacexrResultTests(unittest.TestCase):
    def test_doublet_classes_and_filtered_rows_keep_alignment(self):
        result = dict(barcodes=['s3', 's1', 's0', 's2'], predictions=[
            dict(spot_class='reject'),
            dict(spot_class='doublet_certain', first_type='type1', second_type='type0', weights=[3, 1]),
            dict(spot_class='singlet', first_type='type0'),
            dict(spot_class='doublet_uncertain'),
        ])
        fractions = adapter.read_fractions(result, [f's{i}' for i in range(5)], ['type0', 'type1'], 'doublet')
        np.testing.assert_allclose(fractions.iloc[:2], [[1, 0], [.25, .75]])
        self.assertTrue(fractions.iloc[2:].isna().all().all())

    def test_full_and_multi_normalize_selected_types(self):
        for mode in ['full', 'multi']:
            with self.subTest(mode=mode):
                result = dict(barcodes=['s1', 's0'], predictions=[
                    dict(types=['type2', 'type0'], weights=[2, 6]),
                    dict(types=['type1'], weights=[.8]),
                ])
                fractions = adapter.read_fractions(result, ['s0', 's1', 's2'],
                                                    ['type0', 'type1', 'type2'], mode)
                np.testing.assert_allclose(fractions.iloc[:2], [[0, 1, 0], [.75, 0, .25]])
                self.assertTrue(fractions.iloc[2].isna().all())

    def test_malformed_predictions_are_rejected(self):
        invalid = [
            dict(barcodes=['s0', 's0'], predictions=[]),
            dict(barcodes=['unknown'], predictions=[]),
            dict(barcodes=['s0'], predictions=[]),
        ]
        for types, weights in [(['bad'], [1]), (['type0', 'type0'], [.5, .5]),
                               (['type0'], [-1]), (['type0'], [float('inf')]),
                               (['type0'], [None]), (['type0'], [0]), (['type0'], [1, 2])]:
            invalid.append(dict(barcodes=['s0'], predictions=[dict(types=types, weights=weights)]))
        for result in invalid:
            with self.subTest(result=result), self.assertRaises(ValueError):
                adapter.read_fractions(result, ['s0'], ['type0', 'type1'], 'full')


class SpacexrBridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.target = ad.AnnData(sparse.csr_matrix([[30., 10., 0.], [0., 20., 2.], [2., 2., 2.], [0., 0., 0.]]),
                                obs=pd.DataFrame(index=['001', '002', 'NA', 'zero']),
                                var=pd.DataFrame(index=['g0', 'g1', 'g2']))
        self.reference = ad.AnnData(sparse.csr_matrix([[0., 20., 1.], [0., 30., 1.], [30., 0., 2.], [20., 0., 2.]]),
                                   obs=pd.DataFrame({'cell_type': ['NA', 'NA', 'A/B', 'A/B']},
                                                    index=['001', '002', 'NA', 'with\tspace']),
                                   var=pd.DataFrame(index=['g0', 'g1', 'g2']))

    def job(self, params=None):
        inp, out = self.root / 'input', self.root / 'output'
        inp.mkdir()
        self.target.write_h5ad(inp / 'target.h5ad')
        self.reference.write_h5ad(inp / 'reference.h5ad')
        (inp / 'params.json').write_text(json.dumps(dict(
            params=PLUGIN.resolve_params(params), label_column='cell_type')))
        return sdk.Job(inp, out)

    def fake_r(self, argv, check):
        self.assertTrue(check)
        self.assertEqual(argv[1], '--vanilla')
        bridge = Path(argv[-1])
        spec = json.loads((bridge / 'spec.json').read_text())
        self.assertEqual(spec['reference_labels'], ['type1', 'type1', 'type0', 'type0'])
        self.assertEqual(mmread(bridge / 'reference.mtx').shape, (3, 4))
        result = dict(barcodes=['s1', 's0'], genes=['g0', 'g1'], predictions=[
            dict(types=['type1', 'type0'], weights=[2, 0]),
            dict(types=['type0', 'type1'], weights=[3, 1]),
        ])
        (bridge / 'rctd.json').write_text(json.dumps(result))

    def test_sparse_bridge_restores_exact_labels_and_caps_threads(self):
        self.target.obsm['spatial'] = np.arange(8).reshape(4, 2)
        job = self.job(dict(doublet_mode='full', max_cores=8))
        with patch.object(adapter.subprocess, 'run', side_effect=self.fake_r), \
                patch.dict(os.environ, VHD_DECONV_THREADS='2'):
            adapter.run(job)
        result = DeconvolutionService._read(job.output_dir)
        validate_result(result, self.target.obs_names, ['A/B', 'NA'])
        np.testing.assert_allclose(result.fractions.iloc[:2], [[.75, .25], [0, 1]])
        self.assertTrue(result.fractions.iloc[2:].isna().all().all())
        spec = json.loads((job.output_dir / 'r_inputs' / 'spec.json').read_text())
        self.assertEqual(spec['params']['max_cores'], 2)
        self.assertEqual(spec['coordinates'], self.target.obsm['spatial'].tolist())
        self.assertEqual(json.loads((job.output_dir / 'observation_labels.json').read_text()),
                         dict(zip(['s0', 's1', 's2', 's3'], self.target.obs_names)))
        np.testing.assert_array_equal(mmread(job.output_dir / 'r_inputs' / 'target.mtx').toarray(),
                                      self.target.X.T.toarray())

    def test_fractional_counts_fail_before_starting_r(self):
        self.target.X[0, 0] = .5
        job = self.job()
        with patch.object(adapter.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'raw integer counts'):
            adapter.run(job)
        execute.assert_not_called()

    def test_missing_r_has_setup_guidance(self):
        job = self.job()
        with patch.object(adapter.subprocess, 'run', side_effect=FileNotFoundError), \
                self.assertRaisesRegex(RuntimeError, 'VHD_SPACEXR_RSCRIPT'):
            adapter.run(job)

    def test_service_preserves_spatial_rows_and_publishes_session_cache(self):
        self.target.obsm['spatial'] = np.arange(8).reshape(4, 2)
        directory = self.root / 'session-id' / 'deconvolution'
        service = DeconvolutionService(dict(runner='local', work_dir=str(directory), threads_per_run=1))
        coordinates = []

        def run_plugin(plugin, input_dir, output_dir, progress, log_path):
            job = sdk.Job(input_dir, output_dir)
            target, _ = job.read_inputs()
            coordinates.append(target.obsm['spatial'].copy())
            with patch.object(adapter.subprocess, 'run', side_effect=self.fake_r):
                adapter.run(job)
            log_path.write_text('RCTD bridge test with simulated R output\n')

        with patch.object(service.runner, 'run', side_effect=run_plugin) as execute:
            result = service.run('spacexr', self.target, self.reference, 'cell_type', dict(doublet_mode='full'))
            validate_result(result, self.target.obs_names, ['A/B', 'NA'])
            np.testing.assert_array_equal(coordinates[0], self.target.obsm['spatial'][:3])
            self.assertTrue(Path(result.metadata['log_path']).is_relative_to(directory / 'spacexr'))
            cached = service.run('spacexr', self.target, self.reference, 'cell_type', dict(doublet_mode='full'))
            self.assertEqual(execute.call_count, 1)
            self.assertEqual(cached.metadata['run_id'], result.metadata['run_id'])
            self.target.obsm['spatial'][0, 0] += 1
            changed = service.run('spacexr', self.target, self.reference, 'cell_type', dict(doublet_mode='full'))
            self.assertEqual(execute.call_count, 2)
            self.assertNotEqual(changed.metadata['run_id'], result.metadata['run_id'])

    def test_reference_and_parameter_validation(self):
        service = DeconvolutionService(dict(runner='local'))
        resolved, columns = service.validate_params('spacexr', self.reference, 'cell_type')
        self.assertEqual(columns, ['cell_type'])
        self.assertEqual(resolved['doublet_mode'], 'doublet')
        for params in [dict(UMI_min=20, UMI_max=10), dict(n_max_cells=2), dict(doublet_mode='invalid')]:
            with self.subTest(params=params), self.assertRaises(ValueError):
                service.validate_params('spacexr', self.reference, 'cell_type', params)
        self.reference.obs['cell_type'] = 'only one'
        with self.assertRaisesRegex(ValueError, 'at least two'):
            service.validate_params('spacexr', self.reference, 'cell_type')


if __name__ == '__main__':
    unittest.main()
