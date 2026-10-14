"""Baseline nucleus segmentation: hematoxylin threshold + distance-transform watershed.

Also the reference example for contributors: a complete plugin in ~60 lines.
"""

import numpy as np
from scipy import ndimage as ndi
from skimage.color import rgb2hed
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed

from vhd_seg import count_blocks, iter_blocks, labels_to_polygons, start


def hematoxylin(rgb, sigma):
    channel = rgb2hed(rgb)[..., 0]
    return gaussian(channel, sigma) if sigma > 0 else channel


def segment(tile, threshold, sigma, min_distance, min_area):
    mask = hematoxylin(tile, sigma) > threshold
    mask = remove_small_objects(mask, max(1, int(min_area)))
    distance = ndi.distance_transform_edt(mask)
    components, _ = ndi.label(mask)
    peaks = peak_local_max(distance, min_distance=min_distance, labels=components,
                           exclude_border=False)
    markers = np.zeros(mask.shape, dtype=np.int32)
    markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
    return watershed(-distance, markers, mask=mask)


def main():
    job = start()
    p = job.params

    job.report(0.02, "Reading image")
    image = job.read_image()

    sigma = job.um_to_px(p["smoothing_um"])
    min_distance = max(1, int(round(job.um_to_px(p["min_distance_um"]))))
    min_area = job.um2_to_px2(p["min_area_um2"])

    # One threshold for the whole slide (from a downsampled copy) so tiles agree
    job.report(0.05, "Computing threshold")
    step = max(1, max(image.shape[:2]) // 4000)
    threshold = threshold_otsu(hematoxylin(image[::step, ::step], sigma / step))
    threshold *= p["threshold_factor"]

    overlap = 4 * min_distance + 16
    total = count_blocks(image.shape, p["block_size"])
    polygons = []
    for i, (core, (y0, y1, x0, x1)) in enumerate(iter_blocks(image.shape, p["block_size"], overlap)):
        job.report(0.08 + 0.88 * i / total, f"Segmenting tile {i + 1}/{total}")
        labels = segment(image[y0:y1, x0:x1], threshold, sigma, min_distance, min_area)
        polygons += labels_to_polygons(labels, offset=(y0, x0), core=core, min_area=min_area)

    job.report(0.97, "Writing polygons")
    job.write_polygons(polygons)


if __name__ == "__main__":
    main()
