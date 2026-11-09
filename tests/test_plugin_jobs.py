"""Multi-user safety of plugin runs: timeouts, kills, resource caps, shared jobs, secrets."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from hd_viewer.plugins.jobs import Cancelled, JobTable, JOBS
from hd_viewer.plugins.runners import DockerRunner, ApptainerRunner, LocalRunner, RunnerError
from hd_viewer.analysis.deconvolution import DeconvolutionService, discover as discover_deconv
from hd_viewer.analysis.deconvolution.plugins import ROOT, SDK_DIR
from hd_viewer.startup.data.image.segmentation import SegmentationService
from tests.test_deconvolution import reference
from tests.test_deconvolution_plugins import target

KWARGS = dict(sdk_dir=SDK_DIR, segmentation_root=ROOT)
SECRETS = {'VHD_CIBERSORTX_USERNAME': 'private-user', 'VHD_CIBERSORTX_TOKEN': 'private-token'}


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # A killed child of ours may linger as a zombie until reaped
    try:
        return Path(f'/proc/{pid}/stat').read_text().split()[2] != 'Z'
    except OSError:
        return False


class TimeoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def test_timeout_kills_the_process_group_and_stops_the_run(self):
        pidfile = self.work / 'child.pid'
        # The grandchild stands for what apptainer/R start: it must die too
        script = (f"import subprocess, time; p = subprocess.Popen(['sleep', '60']); "
                  f"open({str(pidfile)!r}, 'w').write(str(p.pid)); print('started', flush=True); time.sleep(60)")
        stopped = []
        runner = LocalRunner({}, **KWARGS)
        log = io.StringIO()
        start = time.monotonic()
        with self.assertRaisesRegex(RunnerError, 'timed out'):
            runner._stream([sys.executable, '-c', script], None, log, timeout=1.5, stop=lambda: stopped.append(1))
        self.assertLess(time.monotonic() - start, 20)
        self.assertEqual(stopped, [1])
        self.assertIn('timed out', log.getvalue())
        pid = int(pidfile.read_text())
        deadline = time.monotonic() + 5
        while alive(pid) and time.monotonic() < deadline:
            time.sleep(.05)
        self.assertFalse(alive(pid))

    def test_timeout_from_config_removes_the_named_container(self):
        calls = self.work / 'calls.txt'
        docker = self.work / 'docker'
        docker.write_text('#!/bin/sh\n'
                          f'echo "$@" >> {calls}\n'
                          'case "$1" in run) exec sleep 60 ;; esac\n'
                          'exit 0\n')
        docker.chmod(0o755)
        runner = DockerRunner({'docker': {'binary': str(docker)}, 'timeout_minutes': 0.02}, **KWARGS)
        plugin = discover_deconv()['nnls']
        with self.assertRaisesRegex(RunnerError, r'timed out after 0\.02 minutes'):
            runner.run(plugin, self.work, self.work, lambda f, m: None, self.work / 'run.log')
        lines = calls.read_text().splitlines()
        name = next(line.split('--name ')[1].split()[0] for line in lines if line.startswith('run '))
        self.assertIn(f'rm -f {name}', lines)
        self.assertTrue(name.startswith('vhd-nnls-'))

    def test_image_inspect_and_pull_have_timeouts(self):
        runner = DockerRunner({'docker': {'image_timeout_minutes': 5}}, **KWARGS)
        plugin = discover_deconv()['nnls']
        with patch('hd_viewer.plugins.runners.subprocess.run',
                   side_effect=subprocess.TimeoutExpired('docker', 60)) as run:
            with self.assertRaisesRegex(RunnerError, 'timed out'):
                runner.prepare(plugin, lambda f, m: None, io.StringIO())
        self.assertIn('timeout', run.call_args.kwargs)
        with patch('hd_viewer.plugins.runners.subprocess.run',
                   return_value=subprocess.CompletedProcess([], 1)), \
                patch.object(runner, '_stream') as stream:
            runner.prepare(plugin, lambda f, m: None, io.StringIO())
        self.assertEqual(stream.call_args.kwargs['timeout'], 300)


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.plugin = discover_deconv()['nnls']
        self.folder = Path(tempfile.gettempdir())

    def test_docker_caps_memory_swap_and_processes(self):
        runner = DockerRunner({'memory': '32g', 'pids_limit': 512}, **KWARGS)
        argv, _ = runner.command(self.plugin, self.folder, self.folder, {}, name='vhd-nnls-x')
        text = ' '.join(map(str, argv))
        for flags in ('--memory 32g', '--memory-swap 32g', '--pids-limit 512', '--name vhd-nnls-x'):
            self.assertIn(flags, text)
        self.assertLess(argv.index('--memory'), argv.index(self.plugin.image))
        argv, _ = DockerRunner({}, **KWARGS).command(self.plugin, self.folder, self.folder, {})
        self.assertFalse({'--memory', '--pids-limit', '--name'} & set(argv))
        self.assertNotIn('--network', argv)

    def test_apptainer_memory_only_when_configured(self):
        argv, _ = ApptainerRunner({'memory': '32g'}, **KWARGS).command(self.plugin, self.folder, self.folder, {})
        self.assertNotIn('--memory', argv)
        argv, _ = ApptainerRunner({'apptainer': {'memory': '16g'}}, **KWARGS).command(
            self.plugin, self.folder, self.folder, {})
        self.assertEqual(argv[argv.index('--memory') + 1], '16g')


class CredentialTests(unittest.TestCase):
    def setUp(self):
        from hd_viewer.analysis.deconvolution.cibersortx import CREDENTIALS
        self.plugin = discover_deconv()['cibersortx']
        self.plugin.command = ['--mixture', 'mixture.tsv']
        self.plugin.secret_args = CREDENTIALS
        self.folder = Path(tempfile.gettempdir())
        patcher = patch.dict(os.environ, SECRETS)
        patcher.start()
        self.addCleanup(patcher.stop)

    def assert_secret_free(self, argv):
        text = ' '.join(map(str, argv))
        self.assertNotIn('private-user', text)
        self.assertNotIn('private-token', text)

    def test_docker_passes_credentials_through_the_environment(self):
        runner = DockerRunner({}, **KWARGS)
        inspected = subprocess.CompletedProcess([], 0, stdout='["/src/CIBERSORTxFractions"]\n', stderr='')
        with patch.object(runner, '_inspect', return_value=inspected):
            argv, env = runner.command(self.plugin, self.folder, self.folder, {})
        self.assert_secret_free(argv)
        self.assertEqual(env, SECRETS)
        for name in SECRETS:
            self.assertEqual(argv[argv.index(name) - 1], '-e')
        self.assertEqual(argv[argv.index('--entrypoint') + 1], '/bin/sh')
        script = argv[argv.index('-c') + 1]
        self.assertEqual(argv[argv.index('-c') + 2:], ['sh', '/src/CIBERSORTxFractions', '--mixture', 'mixture.tsv'])
        # What the shell inside the image runs
        out = subprocess.run(['sh', '-c', script, 'sh', 'echo', '--mixture', 'mixture.tsv'],
                             env={**os.environ, **SECRETS}, capture_output=True, text=True).stdout
        self.assertEqual(out.split(), ['--mixture', 'mixture.tsv', '--username', 'private-user',
                                       '--token', 'private-token'])

    def test_apptainer_passes_credentials_through_apptainerenv(self):
        argv, env = ApptainerRunner({}, **KWARGS).command(self.plugin, self.folder, self.folder, {})
        self.assert_secret_free(argv)
        self.assertEqual(argv[1], 'exec')
        self.assertEqual(env['APPTAINERENV_VHD_CIBERSORTX_TOKEN'], 'private-token')
        self.assertEqual(env['SINGULARITYENV_VHD_CIBERSORTX_USERNAME'], 'private-user')
        self.assertIn('/.singularity.d/runscript', argv[argv.index('-c') + 1])

    def test_cibersortx_bridge_keeps_credentials_out_of_the_command(self):
        from hd_viewer.analysis.deconvolution.cibersortx import run
        from hd_viewer.analysis.deconvolution.api.sdk.vhd_deconv import Job  # noqa: F401
        plugin = discover_deconv()['cibersortx']
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inp, out = root / 'input', root / 'output'
            inp.mkdir(); out.mkdir()
            target()[:2].write_h5ad(inp / 'target.h5ad')
            reference().write_h5ad(inp / 'reference.h5ad')
            (inp / 'params.json').write_text(json.dumps(dict(params=plugin.resolve_params(), label_column='cell_type')))
            runner = DockerRunner({}, **KWARGS)
            with patch.object(runner, 'run', side_effect=RunnerError('stop')) as execute, \
                    self.assertRaises(RunnerError):
                run(plugin, runner, inp, out, lambda f, m: None, root / 'run.log')
        selected = execute.call_args.args[0]
        self.assert_secret_free(selected.command)
        self.assertEqual(set(selected.secret_args.values()), set(SECRETS))


class JobTableTests(unittest.TestCase):
    def test_a_job_survives_while_another_session_awaits_it(self):
        table, kills = JobTable(), []
        entered, release = threading.Barrier(3), threading.Event()

        def wait(owner):
            with table.join('key', owner) as job:
                job.on_cancel(lambda: kills.append(owner)) if owner == 'a' else None
                entered.wait()
                release.wait(5)

        threads = [threading.Thread(target=wait, args=(o,)) for o in ('a', 'b')]
        for thread in threads:
            thread.start()
        entered.wait()
        self.assertEqual(table.cancel_owner('a'), 0)
        self.assertEqual(kills, [])
        self.assertEqual(table.cancel_owner('b'), 1)
        self.assertEqual(kills, ['a'])
        release.set()
        for thread in threads:
            thread.join()
        self.assertEqual(len(table), 0)
        # A closed session cannot start new jobs; others can reuse the key
        with self.assertRaises(Cancelled):
            with table.join('key', 'a'):
                pass
        with table.join('key', 'c') as job:
            self.assertFalse(job.cancelled.is_set())

    def test_withdrawing_stops_the_run_only_when_nobody_else_awaits_it(self):
        table, kills = JobTable(), []
        with table.join('key', 'a') as job, table.join('key', 'b'):
            job.on_cancel(lambda: kills.append(job.reason))
            self.assertFalse(table.withdraw(job, 'a'))          # b still awaits the shared run
            self.assertFalse(job.cancelled.is_set())
            self.assertFalse(table.withdraw(job, 'stranger'))   # not one of its owners
            self.assertTrue(table.withdraw(job, 'b'))
            self.assertEqual(kills, ['cancelled by the user'])
            with self.assertRaisesRegex(Cancelled, 'cancelled by the user'):
                job.check()
        # Withdrawing from a job does not close the session: it can start it again
        with table.join('key', 'a') as again:
            self.assertIsNot(again, job)
            self.assertFalse(again.cancelled.is_set())

    def test_rejoining_after_withdrawing_counts_as_awaiting_again(self):
        table = JobTable()
        with table.join('key', 'a') as job, table.join('key', 'b'):
            table.withdraw(job, 'a')
            with table.join('key', 'a') as same:
                self.assertIs(same, job)
                self.assertEqual(table.cancel_owner('b'), 0)     # a awaits it again
                self.assertFalse(job.cancelled.is_set())

    def test_queued_job_leaves_the_queue_when_cancelled(self):
        table, slots = JobTable(), threading.Semaphore(1)
        slots.acquire()  # another run holds the only slot
        errors = []

        def queued():
            with table.join('key', 'a') as job:
                try:
                    with job.hold(slots, poll=.05):
                        pass
                except Cancelled as exc:
                    errors.append(exc)

        thread = threading.Thread(target=queued)
        thread.start()
        time.sleep(.2)
        table.cancel_owner('a')
        thread.join(5)
        self.assertFalse(thread.is_alive())
        self.assertEqual(len(errors), 1)
        self.assertTrue(slots.acquire(blocking=False) is False)

    def test_closing_the_session_kills_its_running_process(self):
        from hd_viewer.core.session import cleanup_session
        runner, errors = LocalRunner({}, **KWARGS), []
        started = threading.Event()

        def work():
            with JOBS.join(('test', 'running'), 'session-x'):
                try:
                    runner._stream([sys.executable, '-c', 'import time; print("go", flush=True); time.sleep(60)'],
                                   None, io.StringIO(), lambda f, m: started.set())
                except Cancelled as exc:
                    errors.append(exc)

        thread = threading.Thread(target=work)
        thread.start()
        self.assertTrue(started.wait(10))
        start = time.monotonic()
        cleanup_session('session-x')
        thread.join(10)
        self.assertFalse(thread.is_alive())
        self.assertLess(time.monotonic() - start, 10)
        self.assertEqual(len(errors), 1)


class SegmentationSharingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def test_service_is_one_per_process_and_does_not_detect_on_creation(self):
        import hd_viewer.core.pipeline as pipeline
        made = []
        with patch.object(pipeline, '_segmentation_service', None), \
                patch.object(pipeline, 'segmentation_kwargs', return_value=dict(config={'runner': 'local'})), \
                patch('hd_viewer.startup.data.image.segmentation.service.make_runner',
                      side_effect=lambda *a, **k: made.append(1) or LocalRunner({}, **KWARGS)):
            services = []
            threads = [threading.Thread(target=lambda: services.append(pipeline.segmentation_service()))
                       for _ in range(8)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            self.assertEqual(len({id(s) for s in services}), 1)
            service = services[0]
            self.assertIsNone(service.runner_name)
            self.assertEqual(made, [])
            done = threading.Event()
            service.detect_runner(lambda runner: done.set())
            self.assertTrue(done.wait(5))
            self.assertEqual(service.runner_name, 'local')
            self.assertEqual(made, [1])

    def test_identical_concurrent_requests_share_one_run(self):
        service = SegmentationService({'runner': 'local', 'max_parallel': 4})
        calls, barrier = [], threading.Barrier(3)

        class Fake:
            name = 'fake'

            def run(self, plugin, input_dir, output_dir, on_progress, log_path):
                calls.append(Path(input_dir).parent)
                time.sleep(.3)
                log_path.write_text('ok\n')
                ring = [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]
                (output_dir / 'cells.geojson').write_text(json.dumps({'type': 'FeatureCollection', 'features': [
                    {'type': 'Feature', 'geometry': {'type': 'Polygon', 'coordinates': [ring]}, 'properties': {}}]}))

        service._runner = Fake()
        image = self.work / 'image.tif'
        image.write_bytes(b'x')
        results = []

        def request(owner):
            barrier.wait()
            cells, _ = service.run('watershed', image, 'img', self.work, owner=owner)
            results.append(len(cells))

        threads = [threading.Thread(target=request, args=(o,)) for o in ('a', 'b', 'c')]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(results, [1, 1, 1])
        self.assertEqual(len(calls), 1)
        self.assertTrue(calls[0].name.endswith('.running'))
        self.assertFalse(list(self.work.glob('watershed/*.running')))

    def test_failed_run_keeps_only_its_log(self):
        service = SegmentationService({'runner': 'local'})

        class Failing:
            name = 'fake'

            def run(self, plugin, input_dir, output_dir, on_progress, log_path):
                log_path.write_text('boom\n')
                raise RunnerError('exited with code 1')

        service._runner = Failing()
        image = self.work / 'image.tif'
        image.write_bytes(b'x')
        with self.assertRaises(RunnerError) as raised:
            service.run('watershed', image, 'img', self.work)
        self.assertNotIn(str(self.work), str(raised.exception))
        self.assertIn('Run id', str(raised.exception))
        run_dir, = self.work.glob('watershed/*.running')
        self.assertEqual([p.name for p in run_dir.iterdir()], ['run.log'])


class DeconvolutionFailureTests(unittest.TestCase):
    def test_failed_run_removes_inputs_and_hides_server_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            work = Path(folder)
            service = DeconvolutionService(dict(runner='local', work_dir=str(work), threads_per_run=1))

            def fail(plugin, input_dir, output_dir, progress, log_path):
                log_path.write_text('fit failed\n')
                raise RunnerError(f'cannot read {input_dir}/target.h5ad')

            with patch.object(service.runner, 'run', side_effect=fail):
                with self.assertRaises(RunnerError) as raised:
                    service.run('nnls', target(), reference(), 'cell_type', owner='s')
            message = str(raised.exception)
            self.assertNotIn(folder, message)
            self.assertIn('Run id', message)
            run_dir, = (work / 'nnls').glob('*.running-*')
            self.assertEqual([p.name for p in run_dir.iterdir()], ['run.log'])
            self.assertFalse(list(work.rglob('*.h5ad')))
            self.assertEqual(len(JOBS), 0)


if __name__ == '__main__':
    unittest.main()
