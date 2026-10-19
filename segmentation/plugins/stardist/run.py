"""StarDist nucleus segmentation (pretrained H&E model).

Memory: the slide stays uint8. Normalisation percentiles are computed once on a
subsample and applied block by block, and no full-slide label image is built.
Normalising the whole slide as float32 needed ~3 copies of 5+ GB on Visium HD images.
"""

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


# Every SUBSAMPLE-th pixel per axis: percentiles barely change, memory drops 1/SUBSAMPLE²
SUBSAMPLE = 4


def block_normalizer(image, low, high):
    """Per-channel percentile normalisation, applied to each block StarDist reads."""
    from csbdeep.data import Normalizer
    from csbdeep.utils import normalize_mi_ma

    sample = image[::SUBSAMPLE, ::SUBSAMPLE]
    mi = np.percentile(sample, low, axis=(0, 1), keepdims=True)
    ma = np.percentile(sample, high, axis=(0, 1), keepdims=True)

    class SlideNormalizer(Normalizer):
        def before(self, x, axes):
            return normalize_mi_ma(x, mi, ma, dtype=np.float32)

        def after(self, *args, **kwargs):
            raise NotImplementedError

        @property
        def do_after(self):
            return False

    return SlideNormalizer()


def main():
    job = start()
    p = job.params

    job.report(0.02, "Reading image")
    image = job.read_image()

    job.report(0.08, "Computing normalisation")
    normalizer = block_normalizer(image, p["norm_low"], p["norm_high"])

    job.report(0.12, "Loading model")
    model = load_model()

    job.report(0.15, "Detecting nuclei")
    _, polys = model.predict_instances_big(
        image, axes="YXC",
        block_size=p["block_size"], min_overlap=320, context=128,
        prob_thresh=p["prob_thresh"], nms_thresh=p["nms_thresh"],
        n_tiles=(4, 4, 1), normalizer=normalizer,
        labels_out=False,  # only the polygons are used
    )

    job.report(0.95, "Writing polygons")
    # coord: (n, 2, n_rays) with rows (y, x)
    polygons = [list(zip(xs, ys)) for ys, xs in polys["coord"]]
    properties = [{"score": float(s)} for s in polys["prob"]]
    job.write_polygons(polygons, properties)


if __name__ == "__main__":
    main()
