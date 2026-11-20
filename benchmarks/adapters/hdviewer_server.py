"""Serve HDViewer with the benchmark dataset offered as an example (run with .venv).

    python hdviewer_server.py --spec JSON --port 8765

The dataset is registered like the built-in examples (as tests/test_api.py
does). When the spec has cells, the segmentation cache is primed with the
benchmark's cells.geojson under the key the app computes for the "cellpose"
method, so loading cells through the normal UI reuses those outlines instead of
running a segmentation model (segmentation itself is not benchmarked).
"""
import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

DATASET_TITLE = "Benchmark dataset"
SEGMENTATION_PLUGIN = "cellpose"


def example_config(spec):
    from benchmarks.adapters import common
    return dict(expression_path=os.path.abspath(common.matrix_path(spec)),
                positions_path=os.path.abspath(common.positions_path(spec)),
                image_path=os.path.abspath(spec["image_path"]), bins_size=int(spec["bin_size"]))


def prime_segmentation(spec, config):
    """Place cells.geojson where SegmentationService.run looks for a cached result."""
    import hashlib
    from hd_viewer.core.pipeline import image_key, segmentation_service
    from hd_viewer.startup.data.VisiumHD import VisiumHD
    service = segmentation_service()
    plugin = service.plugins[SEGMENTATION_PLUGIN]
    dataset = VisiumHD(img_path=config["image_path"], adata_path=config["expression_path"],
                       positions_path=config["positions_path"], bins_size=config["bins_size"])
    um_per_px = (float(dataset.bins_size) / (dataset.pixels_per_bin * dataset.image_scale)
                 if dataset.pixels_per_bin else None)
    spec_ = {"params": plugin.resolve_params({}), "um_per_px": um_per_px}
    key = hashlib.sha1(json.dumps({"image": image_key(dataset.img_path), "plugin": plugin.name,
                                   "version": plugin.version, **spec_}, sort_keys=True).encode()).hexdigest()[:12]
    out = Path(config["image_path"]).resolve().parent / "segmentation" / plugin.name / key
    out.mkdir(parents=True, exist_ok=True)
    target = out / "cells.geojson"
    source = Path(spec["cells_path"]).resolve()
    if target.resolve() != source:
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(source)
    return plugin.title


def write_map_boxes(spec, config):
    """The shared view boxes (full-res pixels) in map coordinates, for the browser driver.
    Written before the server starts; this preparation is not part of the app's timings."""
    import time
    import numpy as np
    from hd_viewer.core.pipeline import fit_transform
    from hd_viewer.startup.data.VisiumHD import VisiumHD
    started = time.perf_counter()
    dataset = VisiumHD(img_path=config["image_path"], adata_path=config["expression_path"],
                       positions_path=config["positions_path"], bins_size=config["bins_size"])
    transform = fit_transform(dataset)

    def to_map(box):
        x0, x1, y0, y1 = box
        c = transform.points(np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=float))
        return [float(c[:, 0].min()), float(c[:, 0].max()), float(c[:, 1].min()), float(c[:, 1].max())]

    views = spec["viewports"]
    boxes = {name: ([to_map(b) for b in box] if name == "pan" else to_map(box)) for name, box in views.items()}
    Path(spec["boxes_out"]).write_text(json.dumps(dict(boxes=boxes, seconds=time.perf_counter() - started)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--port", type=int, required=True)
    args = parser.parse_args()
    spec = json.loads(args.spec)

    from hd_viewer.startup import datasets
    config = example_config(spec)
    datasets.AVAILABLE_DATASETS.clear()
    datasets.AVAILABLE_DATASETS[DATASET_TITLE] = config
    if spec.get("with_cells"):
        print("segmentation plugin:", prime_segmentation(spec, config), flush=True)
    if spec.get("boxes_out"):
        write_map_boxes(spec, config)

    import uvicorn
    from hd_viewer.api.main import create_app
    uvicorn.run(create_app(), host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
