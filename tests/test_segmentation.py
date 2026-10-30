import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from hd_viewer.startup.data.image.segmentation import SegmentationService, discover, read_polygons
from hd_viewer.startup.data.image.segmentation.plugins import SDK_DIR, Param
from hd_viewer.startup.data.image.segmentation.runners import ApptainerRunner, DockerRunner, RunnerError

sys.path.insert(0, str(SDK_DIR))
from vhd_seg import iter_blocks, labels_to_polygons  # noqa: E402


def synthetic_he(height=600, width=900, seed=0):
    """Pink background with purple disks; returns the image and the disk centres."""
    rng = np.random.default_rng(seed)
    image = np.full((height, width, 3), (235, 200, 225), np.uint8)
    yy, xx = np.mgrid[:height, :width]
    disks = []
    while len(disks) < 60:
        cy, cx, r = rng.integers(15, height - 15), rng.integers(15, width - 15), rng.integers(5, 9)
        if any((cy - a) ** 2 + (cx - b) ** 2 < (r + c + 3) ** 2 for a, b, c in disks):
            continue
        disks.append((cy, cx, r))
        image[(yy - cy) ** 2 + (xx - cx) ** 2 <= r * r] = (70, 40, 140)
    return image, disks


class PluginManifestTests(unittest.TestCase):
    def test_bundled_plugins_are_valid(self):
        plugins = discover()
        self.assertTrue({"stardist", "cellpose", "watershed"} <= set(plugins))
        for plugin in plugins.values():
            plugin.resolve_params()  # defaults pass their own validation

    def test_params_are_coerced_and_validated(self):
        plugin = discover()["cellpose"]
        params = plugin.resolve_params({"diameter_um": "12", "model": "nuclei"})
        self.assertEqual(params["diameter_um"], 12.0)
        with self.assertRaises(ValueError):
            plugin.resolve_params({"diameter_um": 500})
        with self.assertRaises(ValueError):
            plugin.resolve_params({"model": "unknown"})
        with self.assertRaises(ValueError):
            plugin.resolve_params({"not_a_param": 1})

    def test_choice_needs_choices(self):
        with self.assertRaises(ValueError):
            Param(name="model", type="choice")


class TileStitchingTests(unittest.TestCase):
    def test_each_object_is_kept_by_exactly_one_tile(self):
        labels = np.zeros((100, 100), np.int32)
        # Objects straddling the tile boundaries at 50
        labels[45:56, 10:20] = 1
        labels[10:20, 45:56] = 2
        labels[45:56, 45:56] = 3
        labels[80:90, 80:90] = 4

        polygons = []
        for core, (y0, y1, x0, x1) in iter_blocks(labels.shape, 50, overlap=10):
            polygons += labels_to_polygons(labels[y0:y1, x0:x1], offset=(y0, x0), core=core)
        self.assertEqual(len(polygons), 4)

        # Outlines are in full-image coordinates
        object_3 = [p for p in polygons if 40 < np.mean(np.array(p)[:, 0]) < 60
                    and 40 < np.mean(np.array(p)[:, 1]) < 60]
        self.assertEqual(len(object_3), 1)
        xs, ys = np.array(object_3[0]).T
        self.assertEqual((xs.min(), xs.max(), ys.min(), ys.max()), (45, 55, 45, 55))


class LocalRunnerTests(unittest.TestCase):
    def setUp(self):
        self.work = Path(tempfile.mkdtemp())
        self.service = SegmentationService({"runner": "local", "max_parallel": 1})
        self.image, self.disks = synthetic_he()
        self.image_path = self.service.prepare_image(self.image, "synthetic", self.work)

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    def test_watershed_finds_every_nucleus_once_and_caches(self):
        progress = []
        cells, out_dir = self.service.run(
            "watershed", self.image_path, "synthetic", self.work,
            params={"block_size": 512}, um_per_px=0.5,
            on_progress=lambda fraction, message: progress.append(fraction),
        )
        self.assertEqual(len(cells), len(self.disks))
        centroids = cells.geometry.centroid
        for cy, cx, r in self.disks:
            inside = ((centroids.x - cx) ** 2 + (centroids.y - cy) ** 2) < r * r
            self.assertEqual(int(inside.sum()), 1)
        self.assertEqual(progress[-1], 1.0)
        self.assertTrue((out_dir / "run.log").exists())

        cached = []
        again, _ = self.service.run(
            "watershed", self.image_path, "synthetic", self.work,
            params={"block_size": 512}, um_per_px=0.5,
            on_progress=lambda fraction, message: cached.append(message),
        )
        self.assertEqual(cached, ["Loaded cached result"])
        self.assertEqual(len(again), len(cells))

    def test_failure_reports_plugin_output(self):
        plugin = self.service.plugins["watershed"]
        broken = self.work / "broken"
        shutil.copytree(plugin.path, broken)
        (broken / "run.py").write_text("import sys\nprint('model weights not found')\nsys.exit(3)\n")
        plugin.path = broken
        try:
            with self.assertRaises(RunnerError) as ctx:
                self.service.run("watershed", self.image_path, "synthetic", self.work)
        finally:
            self.service.plugins = discover()
        self.assertIn("model weights not found", str(ctx.exception))
        self.assertIn("code 3", str(ctx.exception))

    def test_runs_beyond_max_parallel_wait_their_turn(self):
        import threading
        import time

        events = []

        def run(threshold):
            def on_progress(fraction, message):
                if message.startswith(("Queued", "Starting")) or fraction == 1.0:
                    events.append((time.monotonic(), threshold, message))
            # Different parameters so that neither run is served from the cache
            self.service.run("watershed", self.image_path, "synthetic", self.work,
                             params={"threshold_factor": threshold}, um_per_px=0.5,
                             on_progress=on_progress)

        threads = [threading.Thread(target=run, args=(t,)) for t in (1.0, 1.05)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        starts = sorted((t, run) for t, run, m in events if m.startswith("Starting"))
        ends = {run: t for t, run, m in events if not m.startswith(("Queued", "Starting"))}
        self.assertEqual(sum(m.startswith("Queued") for _, _, m in events), 2)
        self.assertEqual(len(starts), 2)
        first = starts[0][1]
        # With max_parallel = 1 the second run only starts after the first has finished
        self.assertGreaterEqual(starts[1][0], ends[first])

    def test_killed_run_says_why_in_the_error_and_the_log(self):
        plugin = self.service.plugins["watershed"]
        killed = self.work / "killed"
        shutil.copytree(plugin.path, killed)
        # What the kernel does to a process over its memory limit
        (killed / "run.py").write_text(
            "import os, signal\nprint('##progress 0.08 Normalising image', flush=True)\n"
            "os.kill(os.getpid(), signal.SIGKILL)\n")
        plugin.path = killed
        try:
            with self.assertRaises(RunnerError) as ctx:
                self.service.run("watershed", self.image_path, "synthetic", self.work)
        finally:
            self.service.plugins = discover()
        self.assertIn("killed, most likely out of memory", str(ctx.exception))
        log = next(self.work.glob("watershed/*.running/run.log")).read_text()
        self.assertTrue(log.rstrip().endswith("(check the RAM limit of the session or job)"), log)

    def test_thread_budget_reaches_the_plugin(self):
        service = SegmentationService({"runner": "local", "threads_per_run": 3})
        plugin = service.plugins["watershed"]
        probe = self.work / "probe"
        shutil.copytree(plugin.path, probe)
        (probe / "run.py").write_text(
            "import os, sys\n"
            "print(os.environ['OMP_NUM_THREADS'], os.environ['OPENBLAS_NUM_THREADS'])\n"
            "sys.exit(1)\n")
        plugin.path = probe
        with self.assertRaises(RunnerError) as ctx:
            service.run("watershed", self.image_path, "synthetic", self.work)
        self.assertIn("3 3", str(ctx.exception))


class ContainerImageTests(unittest.TestCase):
    """Images are pulled from the registry; building is only a fallback."""

    def setUp(self):
        self.work = Path(tempfile.mkdtemp())
        self.calls = self.work / "calls.txt"
        self.plugin = discover()["watershed"]

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    def fake_docker(self, pull_exit):
        # Images are never present locally; `pull` succeeds or fails as requested
        script = self.work / "docker"
        script.write_text(
            "#!/bin/sh\n"
            f"echo \"$@\" >> {self.calls}\n"
            "case \"$1 $2\" in\n"
            "  'image inspect') exit 1 ;;\n"
            f"  pull*) exit {pull_exit} ;;\n"
            "esac\n"
            "exit 0\n")
        script.chmod(0o755)
        return DockerRunner({"docker": {"binary": str(script)}})

    def prepare(self, runner, plugin=None):
        with open(self.work / "log.txt", "w") as log:
            runner.prepare(plugin or self.plugin, lambda f, m: None, log)
        return self.calls.read_text().splitlines()

    def test_pulls_the_manifest_image(self):
        calls = self.prepare(self.fake_docker(pull_exit=0))
        self.assertEqual(calls[-1], f"pull {self.plugin.image}")
        self.assertFalse(any(c.startswith("build") for c in calls))

    def test_builds_from_dockerfile_when_the_pull_fails(self):
        calls = self.prepare(self.fake_docker(pull_exit=1))
        self.assertTrue(calls[-1].startswith(f"build -t {self.plugin.image}"))

    def test_unpublished_image_without_dockerfile_fails_clearly(self):
        source = self.work / "plugin"
        shutil.copytree(self.plugin.path, source)
        (source / "Dockerfile").unlink()
        self.plugin.path = source
        with self.assertRaises(RunnerError) as ctx:
            self.prepare(self.fake_docker(pull_exit=1))
        self.assertIn(f"could not pull {self.plugin.image}", str(ctx.exception))

    def test_plugins_ignore_packages_from_the_home_directory(self):
        # A pip --user numpy 2 in a mounted ~/.local broke TensorFlow inside the image
        runner = ApptainerRunner({})
        argv, _ = runner.command(self.plugin, self.work, self.work, runner.plugin_env())
        self.assertIn("PYTHONNOUSERSITE=1", argv)
        for plugin in discover().values():
            self.assertIn("ENV PYTHONNOUSERSITE=1", (plugin.path / "Dockerfile").read_text())

    def test_apptainer_runs_the_registry_image(self):
        runner = ApptainerRunner({})
        self.assertEqual(runner.image_uri(self.plugin), f"docker://{self.plugin.image}")
        local = ApptainerRunner({"apptainer": {"images": {"watershed": "/sif/watershed.sif"}}})
        self.assertEqual(local.image_uri(self.plugin), "/sif/watershed.sif")


class ReadPolygonsTests(unittest.TestCase):
    def test_invalid_and_degenerate_outlines(self):
        bowtie = [[0, 0], [10, 10], [10, 0], [0, 10], [0, 0]]
        line = [[0, 0], [5, 0], [10, 0], [0, 0]]
        square = [[20, 20], [30, 20], [30, 30], [20, 30], [20, 20]]
        features = [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
                     "properties": {}} for ring in (bowtie, line, square)]
        with tempfile.NamedTemporaryFile("w", suffix=".geojson", delete=False) as f:
            json.dump({"type": "FeatureCollection", "features": features}, f)
        cells = read_polygons(f.name)
        Path(f.name).unlink()
        self.assertEqual(len(cells), 2)
        self.assertTrue(cells.geometry.is_valid.all())
        self.assertTrue((cells.geometry.geom_type == "Polygon").all())


if __name__ == "__main__":
    unittest.main()
