"""StarDist nucleus segmentation (pretrained H&E model)."""

import os
from pathlib import Path

import numpy as np

from vhd_seg import start

MODEL_NAME = "2D_versatile_he"


def load_model():
    from stardist.models import StarDist2D

    # Images bake the weights in at build time; outside a container they are downloaded
    local = Path(os.environ.get("STARDIST_MODELS", "/opt/models/stardist")) / MODEL_NAME
    if local.exists():
        return StarDist2D(None, name=MODEL_NAME, basedir=str(local.parent))
    return StarDist2D.from_pretrained(MODEL_NAME)


def main():
    job = start()
    p = job.params

    job.report(0.02, "Reading image")
    image = job.read_image()

    job.report(0.08, "Normalising image")
    from csbdeep.utils import normalize
    image = normalize(image.astype(np.float32), p["norm_low"], p["norm_high"], axis=(0, 1))

    job.report(0.12, "Loading model")
    model = load_model()

    job.report(0.15, "Detecting nuclei")
    _, polys = model.predict_instances_big(
        image, axes="YXC",
        block_size=p["block_size"], min_overlap=320, context=128,
        prob_thresh=p["prob_thresh"], nms_thresh=p["nms_thresh"],
        n_tiles=(4, 4, 1),
    )

    job.report(0.95, "Writing polygons")
    # coord: (n, 2, n_rays) with rows (y, x)
    polygons = [list(zip(xs, ys)) for ys, xs in polys["coord"]]
    properties = [{"score": float(s)} for s in polys["prob"]]
    job.write_polygons(polygons, properties)


if __name__ == "__main__":
    main()
