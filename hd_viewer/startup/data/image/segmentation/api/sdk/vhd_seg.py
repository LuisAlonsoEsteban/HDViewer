"""Helpers for segmentation plugins. Copied into every plugin image next to run.py.

The full contract is described in segmentation/README.md. In short, a plugin:

    job = start()                      # reads params.json, parses --input/--output
    image = job.read_image()           # (H, W, 3) RGB array, original orientation
    job.report(0.5, "Halfway")         # optional progress
    job.write_polygons(polygons)       # [[x, y], ...] per cell, pixel coordinates

Only numpy is required to import this module; the image/mask helpers import
tifffile, scipy and opencv lazily.
"""

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

DEFAULT_INPUT = "/data/input"
DEFAULT_OUTPUT = "/data/output"

IMAGE_FILE = "image.tif"
PARAMS_FILE = "params.json"
RESULT_FILE = "cells.geojson"

PROGRESS_PREFIX = "##progress"


class Job:
    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        with open(self.input_dir / PARAMS_FILE) as f:
            spec = json.load(f)
        self.params = spec.get("params", {})
        # Physical scale of the input image; None if the dataset does not provide it
        self.um_per_px = spec.get("um_per_px")

    def read_image(self):
        """The tissue image as (H, W, 3) RGB, row 0 at the top."""
        import tifffile
        return tifffile.imread(self.input_dir / IMAGE_FILE)

    def um_to_px(self, um):
        """Convert a length in µm to pixels (identity when the scale is unknown)."""
        return um / self.um_per_px if self.um_per_px else um

    def um2_to_px2(self, um2):
        return um2 / self.um_per_px ** 2 if self.um_per_px else um2

    def report(self, fraction, message=""):
        """Report progress in [0, 1]; the host shows it on the loading page."""
        fraction = min(1.0, max(0.0, float(fraction)))
        print(f"{PROGRESS_PREFIX} {fraction:.4f} {message}", flush=True)

    def write_polygons(self, polygons, properties=None):
        """
        Write the segmentation result.

        polygons   : iterable of cell outlines, each a sequence of (x, y) vertices in
                     pixel coordinates of the input image (x = column, y = row).
        properties : optional list of dicts (one per polygon), e.g. {"score": 0.93}.
        """
        features = []
        for i, polygon in enumerate(polygons):
            ring = [[float(x), float(y)] for x, y in polygon]
            if len(ring) < 3:
                continue
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            props = properties[i] if properties is not None else {}
            features.append({
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": [ring]},
                "properties": props,
            })

        tmp = self.output_dir / (RESULT_FILE + ".tmp")
        with open(tmp, "w") as f:
            json.dump({"type": "FeatureCollection", "features": features}, f)
        os.replace(tmp, self.output_dir / RESULT_FILE)
        self.report(1.0, f"{len(features)} cells")


def start(argv=None):
    parser = argparse.ArgumentParser(description="Visium HD segmentation plugin")
    parser.add_argument("--input", default=os.environ.get("VHD_SEG_INPUT", DEFAULT_INPUT))
    parser.add_argument("--output", default=os.environ.get("VHD_SEG_OUTPUT", DEFAULT_OUTPUT))
    args = parser.parse_args(argv)
    return Job(args.input, args.output)


def iter_blocks(shape, block_size, overlap):
    """
    Split an image into tiles for methods that cannot process it whole.

    Yields (core, outer) as (y0, y1, x0, x1) boxes: run the method on `outer`
    (core plus `overlap` context) and keep the objects whose centroid lies in
    `core`, so that every cell is kept exactly once.
    """
    height, width = shape[:2]
    for y0 in range(0, height, block_size):
        for x0 in range(0, width, block_size):
            core = (y0, min(y0 + block_size, height), x0, min(x0 + block_size, width))
            outer = (max(0, core[0] - overlap), min(height, core[1] + overlap),
                     max(0, core[2] - overlap), min(width, core[3] + overlap))
            yield core, outer


def count_blocks(shape, block_size):
    return math.ceil(shape[0] / block_size) * math.ceil(shape[1] / block_size)


def labels_to_polygons(labels, offset=(0, 0), core=None, min_area=0):
    """
    Outline every object of an integer label mask (0 = background).

    offset   : (y, x) position of `labels` in the full image.
    core     : optional (y0, y1, x0, x1) box in full-image coordinates; objects
               whose centroid falls outside it are dropped (see iter_blocks).
    min_area : objects smaller than this many pixels are dropped.
    """
    import cv2
    from scipy import ndimage

    oy, ox = offset
    polygons = []
    for label, box in enumerate(ndimage.find_objects(labels), start=1):
        if box is None:
            continue
        mask = (labels[box] == label).astype(np.uint8)
        area = int(mask.sum())
        if area == 0 or area < min_area:
            continue

        top, left = box[0].start + oy, box[1].start + ox
        if core is not None:
            ys, xs = np.nonzero(mask)
            cy, cx = ys.mean() + top, xs.mean() + left
            if not (core[0] <= cy < core[1] and core[2] <= cx < core[3]):
                continue

        # Pad so that objects touching the box edge still get a closed contour
        contours, _ = cv2.findContours(np.pad(mask, 1), cv2.RETR_EXTERNAL,
                                       cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        contour = max(contours, key=cv2.contourArea)[:, 0, :].astype(float)
        if len(contour) < 3:
            continue
        contour[:, 0] += left - 1
        contour[:, 1] += top - 1
        polygons.append(contour.tolist())
    return polygons
