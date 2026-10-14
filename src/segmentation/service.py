"""Run segmentation plugins on a tissue image and load their results."""

import hashlib
import json
import os
import shutil
import threading
from pathlib import Path

import geopandas as gpd
import numpy as np
import yaml
from shapely.geometry import Polygon, shape
from shapely.validation import make_valid

from .plugins import SEGMENTATION_ROOT, discover
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
    def __init__(self, config=None):
        self.config = load_config() if config is None else config
        self.plugins = discover()
        self.runner = make_runner(self.config)
        # Runs beyond max_parallel wait their turn (segmentation is memory hungry)
        self._slots = threading.Semaphore(max(1, int(self.config.get("max_parallel", 2))))
        self._input_lock = threading.Lock()

    def prepare_image(self, image, image_key, work_dir, bgr=False):
        """Write the image once per dataset; every plugin run links to it."""
        path = Path(work_dir) / "_inputs" / image_key / IMAGE_FILE
        with self._input_lock:
            if not path.exists():
                import tifffile
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_suffix(".tmp.tif")
                tifffile.imwrite(tmp, to_rgb(image, bgr=bgr), bigtiff=True, photometric="rgb")
                os.replace(tmp, path)
        return path

    def run(self, name, image_path, image_key, work_dir, params=None, um_per_px=None,
            on_progress=None):
        """
        Segment with plugin `name`. Results are cached per image, plugin version
        and parameters, so re-loading a dataset reuses them.

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

        on_progress(0.0, "Queued: waiting for a free slot")
        with self._slots:
            run_dir = out_dir.with_name(key + ".running")
            shutil.rmtree(run_dir, ignore_errors=True)
            input_dir, output_dir = run_dir / "input", run_dir / "output"
            input_dir.mkdir(parents=True)
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
            except RunnerError as exc:
                raise RunnerError(f"{plugin.title} {exc}\nFull log: {run_dir / LOG_FILE}") from None
            if not (output_dir / RESULT_FILE).exists():
                raise RunnerError(f"{plugin.title} finished without writing {RESULT_FILE}. "
                                  f"Log: {run_dir / LOG_FILE}")

            # Publish atomically: a half-written run is never mistaken for a cached one
            out_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(run_dir / LOG_FILE), out_dir / LOG_FILE)
            for path in output_dir.iterdir():
                if path.name != RESULT_FILE:
                    shutil.move(str(path), out_dir / path.name)
            os.replace(output_dir / RESULT_FILE, result)
            shutil.rmtree(run_dir, ignore_errors=True)

        return read_polygons(result), out_dir
