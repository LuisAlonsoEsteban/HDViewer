"""Plugin contracts, real local NNLS, and viewer persistence/regression tests."""
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from bokeh.models import Button, NumericInput, Select
from bokeh.events import ButtonClick

from hd_viewer.board.operations.layer.deconvolution import DeconvolutionService, DeconvolutionResult, discover, validate_result
from hd_viewer.board.operations.layer.deconvolution.service import matrix_fingerprint
from hd_viewer.board.operations.layer.deconvolution.plugins import ROOT, SDK_DIR
from hd_viewer.board.operations.layer.deconvolution.api.sdk.vhd_deconv import write_table, read_table, Job
from hd_viewer.plugins.manifest import Param
from hd_viewer.plugins.runners import DockerRunner, ApptainerRunner, RunnerError
from tests.test_deconvolution import reference, ImmediateDoc
from tests.test_param_form import make_layer


def target():
    return ad.AnnData(sparse.csr_matrix([[30., 10.], [0., 20.], [0., 0.]]),
                      obs=pd.DataFrame(index=['001', '002', '003']), var=pd.DataFrame(index=['g0', 'g1']))


class PluginServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.service = DeconvolutionService(dict(runner='local', work_dir=str(self.work), threads_per_run=1))

    def test_discovery_and_required_method_parameters(self):
        self.assertEqual(set(self.service.plugins), {'nnls', 'cell2location', 'music', 'spacexr', 'cibersortx'})
        self.assertEqual(self.service.plugins['nnls'].resolve_params()['batch_size'], 256)
        with self.assertRaisesRegex(ValueError, 'required'):
            self.service.validate_params('cell2location', reference(), 'cell_type')
        with self.assertRaisesRegex(ValueError, 'required'):
            self.service.validate_params('music', reference(), 'cell_type')
        ref = reference()
        ref.obs['donor'] = ['d1', 'd1', 'd1']
        with self.assertRaisesRegex(ValueError, 'two reference subjects'):
            self.service.validate_params('music', ref, 'cell_type', {'sample_column': 'donor'})
        ref.obs['donor'] = ['d1', 'd2', 'd2']
        resolved, columns = self.service.validate_params('music', ref, 'cell_type', {'sample_column': 'donor'})
        self.assertEqual(columns, ['cell_type', 'donor'])
        self.assertEqual(resolved['sample_column'], 'donor')
        with self.assertRaisesRegex(ValueError, 'Maximum|maximum'):
            self.service.validate_params('cibersortx', ref, 'cell_type', {'g_min': 1000, 'g_max': 200})

    def test_numeric_and_boolean_validation(self):
        for value in (1.2, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                Param('n', 'int').coerce(value)
        with self.assertRaises(ValueError):
            Param('n', 'float').coerce(float('nan'))
        self.assertEqual(Param('n', 'int', default=256).coerce(None), 256)
        self.assertFalse(Param('flag', 'bool').coerce('false'))
        with self.assertRaises(ValueError):
            Param('flag', 'bool').coerce('invalid')

    def test_real_local_nnls_and_cache_do_not_mutate_inputs(self):
        data, ref = target(), reference()
        counts, annotations = data.X.copy(), ref.obs.copy()
        messages = []
        result = self.service.run('nnls', data, ref, 'cell_type', {'batch_size': 1},
                                  on_progress=lambda f, m: messages.append(m))
        np.testing.assert_allclose(result.fractions.iloc[:2], [[.75, .25], [0, 1]])
        self.assertTrue(result.fractions.iloc[2].isna().all())
        self.assertEqual(result.fractions.index.tolist(), ['001', '002', '003'])
        self.assertTrue(Path(result.metadata['log_path']).exists())
        self.assertIn('Writing target counts for the plugin', messages)
        self.assertIn('Writing single-cell reference for the plugin', messages)
        self.assertTrue(any(message.startswith('Launching NNLS') for message in messages))
        self.assertEqual(messages[-1], 'Deconvolution finished')
        np.testing.assert_array_equal(data.X.toarray(), counts.toarray())
        pd.testing.assert_frame_equal(ref.obs, annotations)
        self.assertFalse(data.obsm)
        with patch.object(self.service.runner, 'run', side_effect=AssertionError('cache miss')):
            cached = self.service.run('nnls', data, ref, 'cell_type', {'batch_size': 1})
        pd.testing.assert_frame_equal(result.fractions, cached.fractions)
        self.assertEqual(result.metadata['run_id'], cached.metadata['run_id'])

    def test_run_storage_is_dataset_specific_without_changing_shared_service(self):
        first_dir = self.work / 'dataset one' / 'deconvolution'
        second_dir = self.work / 'dataset two' / 'deconvolution'
        first = self.service.run('nnls', target(), reference(), 'cell_type', work_dir=first_dir)
        first_log = Path(first.metadata['log_path'])
        self.assertTrue(first_log.is_relative_to(first_dir))
        self.assertTrue(first_log.exists())
        self.assertTrue((first_log.parent / 'fractions.json').exists())
        self.assertTrue((first_log.parent / 'result.json').exists())
        self.assertTrue(list((first_dir / '_locks').glob('*.lock')))
        with patch.object(self.service.runner, 'run', side_effect=AssertionError('cache missed')):
            cached = self.service.run('nnls', target(), reference(), 'cell_type', work_dir=first_dir)
        self.assertEqual(cached.metadata['log_path'], str(first_log))
        second = self.service.run('nnls', target(), reference(), 'cell_type', work_dir=second_dir)
        self.assertTrue(Path(second.metadata['log_path']).is_relative_to(second_dir))
        self.assertEqual(first.metadata['run_id'], second.metadata['run_id'])
        self.assertEqual(self.service.work_dir, self.work)
        self.assertFalse((self.work / 'nnls').exists())

    def test_moved_cache_is_reused_and_log_path_follows_new_directory(self):
        first = self.service.run('nnls', target(), reference(), 'cell_type')
        destination = self.work / 'dataset' / 'deconvolution'
        destination.mkdir(parents=True)
        shutil.move(str(self.work / 'nnls'), destination / 'nnls')
        messages = []
        with patch.object(self.service.runner, 'run', side_effect=AssertionError('moved cache was not reused')):
            result = self.service.run('nnls', target(), reference(), 'cell_type', work_dir=destination,
                                      on_progress=lambda fraction, message: messages.append(message))
        self.assertEqual(result.metadata['run_id'], first.metadata['run_id'])
        self.assertTrue(Path(result.metadata['log_path']).is_relative_to(destination))
        self.assertTrue(Path(result.metadata['log_path']).exists())
        self.assertIn('Loaded cached result', messages)

    def test_fingerprint_dense_sparse_and_expression_identity(self):
        data = target()
        dense = data.copy()
        dense.X = data.X.toarray()
        self.assertEqual(matrix_fingerprint(data), matrix_fingerprint(dense))
        dense.X[0, 0] += 1
        self.assertNotEqual(matrix_fingerprint(data), matrix_fingerprint(dense))
        dense = data.copy()
        dense.obs_names = ['003', '002', '001']
        self.assertNotEqual(matrix_fingerprint(data), matrix_fingerprint(dense))

    def test_cache_invalidates_counts_labels_params_and_version(self):
        data, ref = target(), reference()
        first = self.service.run('nnls', data, ref, 'cell_type')
        ref.obs['unrelated'] = ['x', 'y', 'z']
        data.obs['unrelated'] = 1
        self.assertEqual(first.metadata['run_id'], self.service.run('nnls', data, ref, 'cell_type').metadata['run_id'])
        variants = []
        changed = data.copy()
        changed.X[0, 0] += 1
        variants.append(self.service.run('nnls', changed, ref, 'cell_type'))
        changed_ref = ref.copy()
        changed_ref.obs['cell_type'] = ['AA', 'AA', 'B']
        variants.append(self.service.run('nnls', data, changed_ref, 'cell_type'))
        variants.append(self.service.run('nnls', data, ref, 'cell_type', {'batch_size': 1}))
        self.service.plugins['nnls'].version = 'next'
        variants.append(self.service.run('nnls', data, ref, 'cell_type'))
        self.assertEqual(len({first.metadata['run_id'], *[r.metadata['run_id'] for r in variants]}), 5)

    def test_all_zero_target_skips_runner(self):
        data = target()
        data.X = sparse.csr_matrix(data.shape)
        with patch.object(self.service.runner, 'run', side_effect=AssertionError('empty input launched')):
            result = self.service.run('nnls', data, reference(), 'cell_type')
        self.assertTrue(result.fractions.isna().all().all())

    def test_failed_or_invalid_plugin_result_is_never_published(self):
        def invalid(plugin, input_dir, output_dir, progress, log_path):
            log_path.write_text('Invalid adapter output\n')
            frame = pd.DataFrame([[.8, .8], [.5, .5]], index=['001', '002'], columns=['A', 'B'])
            write_table(frame, output_dir / 'fractions.json')
            (output_dir / 'result.json').write_text(json.dumps(dict(genes=['g0', 'g1'])))
        with patch.object(self.service.runner, 'run', side_effect=invalid):
            with self.assertRaisesRegex(RunnerError, 'sum to one'):
                self.service.run('nnls', target(), reference(), 'cell_type')
        self.assertFalse(list((self.work / 'nnls').glob('*/result.json')))
        self.assertTrue(list((self.work / 'nnls').glob('*.running-*/run.log')))
        # The same specification can be retried successfully.
        result = self.service.run('nnls', target(), reference(), 'cell_type')
        self.assertEqual(result.metadata['method'], 'nnls')

    def test_result_validation_rejects_alignment_duplicates_nan_and_abundance(self):
        good = pd.DataFrame([[.75, .25], [0, 1]], index=['001', '002'], columns=['A', 'B'])
        for frame in (good.iloc[::-1], good.set_axis(['A', 'A'], axis=1),
                      good * -1, good * .5, good.assign(A=[np.inf, 0]), good.assign(A=[np.nan, 0])):
            with self.assertRaises(ValueError):
                validate_result(DeconvolutionResult(frame, ['g0']), good.index)
        abundance = good * 10
        validate_result(DeconvolutionResult(good, ['g0'], abundance, {'abundance_units': 'cells'}), good.index)
        abundance.iloc[0, 0] = 0
        with self.assertRaisesRegex(ValueError, 'normalized abundance'):
            validate_result(DeconvolutionResult(good, ['g0'], abundance, {'abundance_units': 'cells'}), good.index)

    def test_reference_artifacts_reused_across_different_target_layers(self):
        calls = []
        def fit(plugin, input_dir, output_dir, progress, log_path):
            job = Job(input_dir, output_dir)
            data, _ = job.read_inputs()
            calls.append((job.reference_cache / 'signatures.json').exists())
            job.reference_output.joinpath('signatures.json').write_text('{}')
            job.write_result(pd.DataFrame([[.5, .5]] * data.n_obs, index=data.obs_names, columns=['A', 'B']), ['g0', 'g1'])
            log_path.write_text('fit\n')
        data = target()
        directory = self.work / 'dataset' / 'deconvolution'
        with patch.object(self.service.runner, 'run', side_effect=fit):
            self.service.run('nnls', data, reference(), 'cell_type', work_dir=directory)
            data.X[0, 0] += 1
            self.service.run('nnls', data, reference(), 'cell_type', work_dir=directory)
        self.assertEqual(calls, [False, True])
        self.assertTrue((directory / '_references' / 'nnls').exists())
        self.assertFalse((self.work / '_references').exists())

    def test_independent_runner_environment_and_cibersortx_requirements(self):
        with patch.dict(os.environ, {'VHD_SEG_RUNNER': 'bad-value'}, clear=False):
            service = DeconvolutionService(dict(runner='local', work_dir=str(self.work)))
        self.assertEqual(service.runner.name, 'local')
        self.assertIn('Docker or Apptainer', service.unavailable_reason('cibersortx'))
        with patch.object(service.runner, 'name', 'docker'), patch.dict(os.environ, {}, clear=True):
            self.assertIn('VHD_CIBERSORTX_TOKEN', service.unavailable_reason('cibersortx'))

    def test_reference_lock_wait_is_reported_and_released(self):
        import fcntl
        import threading
        from hd_viewer.board.operations.layer.deconvolution.service import _file_lock
        path = self.work / 'reference.lock'
        waiting, acquired = threading.Event(), threading.Event()
        with path.open('a') as owner:
            fcntl.flock(owner, fcntl.LOCK_EX)
            def worker():
                with _file_lock(path, on_wait=waiting.set):
                    acquired.set()
            thread = threading.Thread(target=worker)
            thread.start()
            try:
                self.assertTrue(waiting.wait(5))
                self.assertFalse(acquired.is_set())
            finally:
                fcntl.flock(owner, fcntl.LOCK_UN)
                thread.join(timeout=5)
        self.assertFalse(thread.is_alive())
        self.assertTrue(acquired.is_set())

    def test_concurrent_identical_runs_fit_once(self):
        import threading
        calls, results, errors = [], [], []
        start = threading.Barrier(2)
        def fit(plugin, input_dir, output_dir, progress, log_path):
            job = Job(input_dir, output_dir)
            data, _ = job.read_inputs()
            calls.append(1)
            job.write_result(pd.DataFrame([[.5, .5]] * data.n_obs, index=data.obs_names,
                                          columns=['A', 'B']), ['g0', 'g1'])
            log_path.write_text('one fit\n')
        def worker():
            try:
                start.wait()
                results.append(self.service.run('nnls', target(), reference(), 'cell_type'))
            except Exception as exc:
                errors.append(exc)
        with patch.object(self.service.runner, 'run', side_effect=fit):
            threads = [threading.Thread(target=worker) for _ in range(2)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=10)
                self.assertFalse(thread.is_alive())
        self.assertFalse(errors, errors)
        self.assertEqual(len(calls), 1)
        self.assertEqual(results[0].metadata['run_id'], results[1].metadata['run_id'])

    def test_observations_without_shared_expression_remain_unassigned(self):
        data = target()
        extended = ad.AnnData(sparse.hstack([data.X, sparse.csr_matrix([[0], [0], [10]])], format='csr'),
                              obs=data.obs.copy(), var=pd.DataFrame(index=['g0', 'g1', 'extra']))
        result = self.service.run('nnls', extended, reference(), 'cell_type')
        self.assertTrue(result.fractions.iloc[2].isna().all())

    def test_nnls_reports_through_original_progress_utilities(self):
        from hd_viewer.startup.utils.progress import reporting
        from hd_viewer.board.operations.layer.deconvolution import deconvolve
        messages = []
        with reporting(messages.append):
            deconvolve(target(), reference(), 'cell_type')
        self.assertTrue(any('Building' in message for message in messages))
        self.assertTrue(any('Deconvolved' in message for message in messages))

    def test_no_frontend_imports_in_standalone_package(self):
        code = '''import sys
import hd_viewer.board.operations.layer.deconvolution
assert not any(m.startswith(('bokeh', 'hd_viewer.board.layout', 'hd_viewer.board.components', 'hd_viewer.startup')) for m in sys.modules)
'''
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


class ContainerCommandTests(unittest.TestCase):
    def test_unstructured_container_startup_output_is_visible(self):
        runner = DockerRunner({}, sdk_dir=SDK_DIR, segmentation_root=ROOT)
        log, messages = io.StringIO(), []
        runner._stream([sys.executable, '-c', "print('Downloading image layers')"], None, log,
                       lambda fraction, message: messages.append(message))
        self.assertEqual(messages, ['Downloading image layers'])
        self.assertIn('Downloading image layers', log.getvalue())

    def test_deconvolution_images_use_ghcr_for_docker_and_apptainer(self):
        with tempfile.TemporaryDirectory() as folder:
            kwargs = dict(sdk_dir=SDK_DIR, segmentation_root=ROOT)
            apptainer = ApptainerRunner({'apptainer': {'images_dir': folder}}, **kwargs)
            docker = DockerRunner({}, **kwargs)
            for name in ('nnls', 'cell2location', 'music', 'spacexr'):
                plugin = discover()[name]
                expected = f'ghcr.io/luisalonsoesteban/visium-deconv-{name}:0.1.0'
                self.assertEqual(plugin.image, expected)
                self.assertEqual(apptainer.image_uri(plugin), 'docker://' + expected)
                command, _ = docker.command(plugin, Path(folder), Path(folder), docker.plugin_env())
                self.assertEqual(command[-1], expected)

    def test_official_cibersortx_mounts_commands_and_redaction(self):
        plugin = discover()['cibersortx']
        plugin.command = ['--username', 'private-user', '--token', 'private-token']
        kwargs = dict(sdk_dir=SDK_DIR, segmentation_root=ROOT)
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {
                'VHD_CIBERSORTX_USERNAME': 'private-user', 'VHD_CIBERSORTX_TOKEN': 'private-token'}):
            for runner in (DockerRunner({}, **kwargs), ApptainerRunner({}, **kwargs)):
                argv, _ = runner.command(plugin, Path(folder), Path(folder), runner.plugin_env())
                self.assertTrue(any('/src/data:ro' in str(a) for a in argv))
                self.assertTrue(any('/src/outdir' in str(a) for a in argv))
                self.assertEqual(argv[-4:], plugin.command)
                redacted = runner.redact(' '.join(str(a) for a in argv))
                self.assertNotIn('private-token', redacted)
                self.assertNotIn('private-user', redacted)
            runner = DockerRunner({}, **kwargs)
            log = io.StringIO()
            with self.assertRaises(RunnerError) as raised:
                runner._stream([sys.executable, '-c', "print('private-user private-token'); raise SystemExit(3)"], None, log)
            self.assertNotIn('private-token', log.getvalue())
            self.assertNotIn('private-token', str(raised.exception))

    def test_cibersortx_signature_workflow_and_cache(self):
        from hd_viewer.board.operations.layer.deconvolution.cibersortx import run
        plugin = discover()['cibersortx']
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {
                'VHD_CIBERSORTX_USERNAME': 'private-user', 'VHD_CIBERSORTX_TOKEN': 'private-token'}):
            root = Path(folder)
            inp, out = root / 'input', root / 'output'
            inp.mkdir(); out.mkdir()
            target()[:2].write_h5ad(inp / 'target.h5ad')
            reference().write_h5ad(inp / 'reference.h5ad')
            (inp / 'params.json').write_text(json.dumps(dict(params=plugin.resolve_params(), label_column='cell_type')))
            runner = DockerRunner({}, sdk_dir=SDK_DIR, segmentation_root=ROOT)
            def official(selected, input_dir, output_dir, progress, log_path):
                self.assertIn('--single_cell', selected.command)
                self.assertIn('--refsample', selected.command)
                pd.DataFrame([[.75, .25], [0, 1]], index=['s0', 's1'], columns=['type0', 'type1']).to_csv(out / 'CIBERSORTx_Results.txt', sep='\t')
                pd.DataFrame([[1, 0], [0, 1]], index=['g0', 'g1'], columns=['type0', 'type1']).to_csv(out / 'CIBERSORTx_signature.txt', sep='\t')
                log_path.write_text('redacted\n')
            with patch.object(runner, 'run', side_effect=official):
                run(plugin, runner, inp, out, lambda f, m: None, root / 'run.log')
            self.assertTrue((out / 'reference_artifacts' / 'signature.tsv').exists())
            fractions = read_table(out / 'fractions.json')
            self.assertEqual(fractions.index.tolist(), ['001', '002'])
            self.assertEqual(fractions.columns.tolist(), ['A', 'B'])
            # Second target uses the generated signature rather than rebuilding it.
            shutil.copytree(out / 'reference_artifacts', inp / 'reference_artifacts')
            with patch.object(runner, 'run') as execute:
                run(plugin, runner, inp, out, lambda f, m: None, root / 'run.log')
            selected = execute.call_args.args[0]
            self.assertIn('--sigmatrix', selected.command)
            self.assertNotIn('--single_cell', selected.command)


class ViewerPluginTests(unittest.TestCase):
    def test_method_controls_and_failure_restore(self):
        layer, _ = make_layer()
        manager = layer.deconvolution_manager
        manager.reference = reference()
        controls = manager.get_widgets(ImmediateDoc())
        method = next(w for w in controls if isinstance(w, Select) and w.title == 'Deconvolution method')
        run = next(w for w in controls if isinstance(w, Button))
        self.assertEqual(method.value, 'nnls')
        self.assertFalse(run.disabled)
        method.value = 'cell2location'
        self.assertTrue(run.disabled)
        expected = next(w for w in controls[4].children if isinstance(w, NumericInput) and 'Expected cells' in w.title)
        expected.value = 5
        self.assertFalse(run.disabled)
        method.value = 'nnls'
        with patch.object(manager.service, 'run', side_effect=RunnerError('synthetic failure')):
            run._trigger_event(ButtonClick(run))
            import time
            deadline = time.monotonic() + 5
            while run.disabled and time.monotonic() < deadline:
                time.sleep(.01)
        self.assertFalse(run.disabled)
        self.assertFalse(method.disabled)
        self.assertIn('synthetic failure', controls[controls.index(run) + 1].text)
        self.assertFalse(manager.history)

    def test_history_abundance_export_and_saved_readout(self):
        layer, _ = make_layer()
        manager = layer.deconvolution_manager
        manager.reference = reference()
        index = layer.sample.adata.obs_names
        fractions = pd.DataFrame([[.25, .75]] * len(index), index=index, columns=['A', 'B'])
        manager.apply_result(DeconvolutionResult(fractions, ['g0', 'g1'], fractions * 4,
                            dict(method='cell2location', params={'batch_column': None}, abundance_units='expected cells')), 'cell_type')
        manager.apply_result(DeconvolutionResult(fractions, ['g0'], metadata=dict(method='nnls')), 'cell_type')
        self.assertEqual(len(manager.history), 2)
        self.assertNotIn('deconvolution_abundance', layer.sample.adata.obsm)
        self.assertIn('deconvolution_1_abundance', layer.sample.adata.obsm)
        self.assertNotIn('deconv_1_abundance:0', layer.sample.adata.uns['mean_obs_columns'])
        controls = manager.get_widgets(ImmediateDoc())
        selects = {w.title: w for w in controls if isinstance(w, Select)}
        selects['Saved deconvolution run'].value = '1'
        selects['Readout'].value = 'abundance'
        self.assertEqual(selects['Display cell type'].value, 'deconv_1_abundance:0')
        self.assertTrue(layer.sample.is_continous('deconv_1_abundance:0'))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'export.h5ad'
            layer.sample.adata.write_h5ad(path)
            restored = ad.read_h5ad(path)
        self.assertEqual(restored.uns['deconvolution_runs']['1']['method'], 'cell2location')
        self.assertIsNone(json.loads(restored.uns['deconvolution_runs']['1']['params_json'])['batch_column'])
        pd.testing.assert_frame_equal(restored.obsm['deconvolution_1_abundance'], fractions * 4)
        # Restored AnnData history can be browsed without the original in-memory manager.
        layer.sample.adata = restored
        manager.current_run = ''
        self.assertEqual(len(manager.get_widgets(ImmediateDoc())), len(controls))


class ExternalAdapterSmokeTests(unittest.TestCase):
    """Use configured isolated environments; skip missing external dependencies."""
    def smoke(self, name, params):
        python = os.environ.get('VHD_TEST_' + name.upper() + '_PYTHON')
        if not python:
            self.skipTest(f'Set VHD_TEST_{name.upper()}_PYTHON to an isolated plugin interpreter')
        rng = np.random.default_rng(3)
        genes = [f'g{i}' for i in range(30)]
        values = rng.poisson(5, size=(80, 30)).astype(float)
        values[:40, :15] += 20
        values[40:, 15:] += 20
        ref = ad.AnnData(sparse.csr_matrix(values), obs=pd.DataFrame({
            'cell_type': ['A'] * 40 + ['B'] * 40, 'donor': ['d1', 'd2'] * 40},
            index=[f'r{i}' for i in range(80)]), var=pd.DataFrame(index=genes))
        data = ad.AnnData(sparse.csr_matrix(rng.poisson(10, size=(3, 30))),
                          obs=pd.DataFrame(index=[f's{i}' for i in range(3)]), var=ref.var.copy())
        with tempfile.TemporaryDirectory() as folder:
            service = DeconvolutionService(dict(runner='local', work_dir=folder,
                                                local={'python': {name: python}}, threads_per_run=1))
            result = service.run(name, data, ref, 'cell_type', params)
            validate_result(result, data.obs_names, ['A', 'B'])

    def test_cell2location(self):
        self.smoke('cell2location', dict(N_cells_per_location=5, reference_epochs=2,
                                         spatial_epochs=2, posterior_samples=2, batch_size=40))

    def test_music(self):
        self.smoke('music', dict(sample_column='donor', batch_size=2))

    def test_spacexr(self):
        for mode in ('doublet', 'full', 'multi'):
            with self.subTest(mode=mode):
                self.smoke('spacexr', dict(doublet_mode=mode, MAX_MULTI_TYPES=2, UMI_min_sigma=100))

    def test_cibersortx_credentials_and_real_container(self):
        runner = os.environ.get('VHD_TEST_CIBERSORTX_RUNNER')
        if runner not in ('docker', 'apptainer'):
            self.skipTest('Set VHD_TEST_CIBERSORTX_RUNNER=docker or apptainer for the credentialed smoke test')
        if not all(os.environ.get(k) for k in ('VHD_CIBERSORTX_USERNAME', 'VHD_CIBERSORTX_TOKEN')):
            self.skipTest('CIBERSORTx server credentials are unavailable')
        rng = np.random.default_rng(1)
        counts = rng.poisson(5, size=(120, 100)).astype(float)
        counts[:60, :50] += 50
        counts[60:, 50:] += 50
        ref = ad.AnnData(sparse.csr_matrix(counts), obs=pd.DataFrame({
            'cell_type': ['A'] * 60 + ['B'] * 60}, index=[f'r{i}' for i in range(120)]),
            var=pd.DataFrame(index=[f'g{i}' for i in range(100)]))
        data = ad.AnnData(sparse.csr_matrix(np.stack([counts[:60].mean(axis=0), counts[60:].mean(axis=0)])),
                          obs=pd.DataFrame(index=['s0', 's1']), var=ref.var.copy())
        with tempfile.TemporaryDirectory() as folder:
            service = DeconvolutionService(dict(runner=runner, work_dir=folder, threads_per_run=1))
            params = dict(g_min=5, g_max=20, replicates=3)
            result = service.run('cibersortx', data, ref, 'cell_type', params)
            validate_result(result, data.obs_names, ['A', 'B'])
            # A different layer must reuse the signature through the official container.
            data.X[0, 0] += 1
            again = service.run('cibersortx', data, ref, 'cell_type', params)
            validate_result(again, data.obs_names, ['A', 'B'])
            self.assertEqual(result.metadata['reference_id'], again.metadata['reference_id'])
