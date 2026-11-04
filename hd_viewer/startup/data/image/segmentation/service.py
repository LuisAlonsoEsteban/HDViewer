"""Run segmentation plugins on a tissue image and load their results."""


import hashlib
import json
import os
import shutil
import tempfile
import threading
from pathlib import Path

import geopandas as gpd
import numpy as np
import yaml
from shapely.geometry import Polygon, shape
from shapely.validation import make_valid

from .....plugins.jobs import JOBS
from .plugins import SEGMENTATION_ROOT, PLUGINS_DIR, SDK_DIR, discover
from .runners import RunnerError, make_runner

CONFIG_PATH = SEGMENTATION_ROOT / "config.yaml"
IMAGE_FILE = "image.tif"
PARAMS_FILE = "params.json"
RESULT_FILE = "cells.geojson"
LOG_FILE = "run.log"


def load_config(path=CONFIG_PATH):
    if not Path(path).exists():
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


def to_rgb(image, bgr=False):
    """(H, W, 3) RGB view of an image read by tifffile (RGB) or OpenCV (BGR)."""
    image = np.asarray(image)
    if image.ndim == 2:
        image = np.repeat(image[..., None], 3, axis=2)
    if image.shape[2] > 3:
        image = image[..., :3]
    if bgr:
        image = image[..., ::-1]
    return np.ascontiguousarray(image)


def read_polygons(path):
    """cells.geojson -> GeoDataFrame of valid polygons in input-image pixels."""
    with open(path) as f:
        features = json.load(f).get("features", [])
    geometries, properties = [], []
    for feature in features:
        geometry = shape(feature["geometry"])
        if not geometry.is_valid:
            geometry = make_valid(geometry)
            # A self-intersecting outline may become several parts: keep the largest
            if geometry.geom_type != "Polygon":
                parts = [g for g in getattr(geometry, "geoms", []) if g.geom_type == "Polygon"]
                if not parts:
                    continue
                geometry = max(parts, key=lambda g: g.area)
        if geometry.is_empty or geometry.area == 0:
            continue
        geometries.append(geometry)
        properties.append(feature.get("properties") or {})
    gdf = gpd.GeoDataFrame(properties, geometry=geometries)
    return gdf


class SegmentationService:
    def __init__(self, config=None, *, config_path=CONFIG_PATH, plugins_dir=PLUGINS_DIR,
                 sdk_dir=SDK_DIR, segmentation_root=SEGMENTATION_ROOT):
        self.config = load_config(config_path) if config is None else config
        self.plugins = discover(plugins_dir)
        self._runner_args = dict(sdk_dir=sdk_dir, segmentation_root=segmentation_root)
        self._runner = None
        self._runner_lock = threading.Lock()
        # Runs beyond max_parallel wait their turn (segmentation is memory hungry)
        self._slots = threading.Semaphore(max(1, int(self.config.get("max_parallel", 2))))

    @property
    def runner(self):
        """Detected on first use: `auto` probes the Docker daemon, which may take seconds."""
        with self._runner_lock:
            if self._runner is None:
                self._runner = make_runner(self.config, **self._runner_args)
            return self._runner

    @property
    def runner_name(self):
        """The runner's name, or None while it is not detected yet (never blocks)."""
        return self._runner.name if self._runner is not None else None

    def detect_runner(self, callback):
        """Detect the runner in a background thread; callback(runner or exception)."""
        def detect():
            try:
                result = self.runner
            except Exception as exc:
                result = exc
            callback(result)
        threading.Thread(target=detect, daemon=True).start()

    def prepare_image(self, image, image_key, work_dir, bgr=False):
        """Write the image once per dataset; every plugin run links to it."""
        path = Path(work_dir) / "_inputs" / image_key / IMAGE_FILE
        # Datasets (e.g. examples) are shared by sessions: one writer per image
        with JOBS.join(("segmentation-input", str(path))) as job, job.lock:
            if not path.exists():
                import tifffile
                path.parent.mkdir(parents=True, exist_ok=True)
                fd, tmp = tempfile.mkstemp(prefix=".image-", suffix=".tif", dir=path.parent)
                os.close(fd)
                try:
                    tifffile.imwrite(tmp, to_rgb(image, bgr=bgr), bigtiff=True, photometric="rgb")
                    os.replace(tmp, path)
                finally:
                    if os.path.exists(tmp):
                        os.remove(tmp)
        return path

    def run(self, name, image_path, image_key, work_dir, params=None, um_per_px=None,
            on_progress=None, owner=None):
        """
        Segment with plugin `name`. Results are cached per image, plugin version
        and parameters, so re-loading a dataset reuses them. Concurrent identical
        requests share one run, which is cancelled once no `owner` (session) awaits it.

        Returns a GeoDataFrame of cells in input-image pixel coordinates and the
        output directory (which also holds run.log).
        """
        on_progress = on_progress or (lambda fraction, message: None)
        plugin = self.plugins[name]
        resolved = plugin.resolve_params(params)
        spec = {"params": resolved, "um_per_px": um_per_px}

        key = hashlib.sha1(json.dumps(
            {"image": image_key, "plugin": plugin.name, "version": plugin.version, **spec},
            sort_keys=True,
        ).encode()).hexdigest()[:12]
        out_dir = Path(work_dir) / plugin.name / key
        result = out_dir / RESULT_FILE
        if result.exists():
            on_progress(1.0, "Loaded cached result")
            return read_polygons(result), out_dir

        with JOBS.join(("segmentation", str(out_dir.resolve())), owner) as job:
            if job.lock.locked():
                on_progress(0.0, "Waiting for the same segmentation started by another session")
            with job.lock:
                if result.exists():
                    on_progress(1.0, "Loaded cached result")
                    return read_polygons(result), out_dir
                job.check()
                on_progress(0.0, "Queued: waiting for a free slot")
                with job.hold(self._slots):
                    self._run(plugin, key, spec, image_path, out_dir, on_progress)
        return read_polygons(result), out_dir

    def _run(self, plugin, key, spec, image_path, out_dir, on_progress):
        result = out_dir / RESULT_FILE
        out_dir.parent.mkdir(parents=True, exist_ok=True)
        # Unique: another viewer process may run the same key
        run_dir = Path(tempfile.mkdtemp(prefix=key + "-", suffix=".running", dir=out_dir.parent))
        input_dir, output_dir = run_dir / "input", run_dir / "output"
        input_dir.mkdir()
        output_dir.mkdir()

        with open(input_dir / PARAMS_FILE, "w") as f:
            json.dump(spec, f, indent=2)
        try:
            os.link(image_path, input_dir / IMAGE_FILE)
        except OSError:
            shutil.copy(image_path, input_dir / IMAGE_FILE)

        on_progress(0.0, f"Starting ({self.runner.name})")
        try:
            self.runner.run(plugin, input_dir, output_dir, on_progress, run_dir / LOG_FILE)
            if not (output_dir / RESULT_FILE).exists():
                raise RunnerError(f"finished without writing {RESULT_FILE}.")
        except BaseException as exc:
            # Failed runs keep only their log (no server paths in the browser)
            shutil.rmtree(input_dir, ignore_errors=True)
            shutil.rmtree(output_dir, ignore_errors=True)
            print(f"[segmentation] {plugin.title} failed, log: {run_dir / LOG_FILE}")
            if isinstance(exc, RunnerError):
                raise RunnerError(f"{plugin.title} {exc}\nRun id: {run_dir.name}") from None
            raise

        # Publish atomically: a half-written run is never mistaken for a cached one
        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.move(str(run_dir / LOG_FILE), out_dir / LOG_FILE)
        for path in output_dir.iterdir():
            if path.name != RESULT_FILE:
                shutil.move(str(path), out_dir / path.name)
        os.replace(output_dir / RESULT_FILE, result)
        shutil.rmtree(run_dir, ignore_errors=True)
