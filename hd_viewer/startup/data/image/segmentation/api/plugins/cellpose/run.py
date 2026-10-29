"""Cellpose segmentation, run tile by tile so that whole-slide images fit in memory."""

import os

from vhd_seg import count_blocks, iter_blocks, labels_to_polygons, start


def main():
    job = start()
    p = job.params

    job.report(0.02, "Reading image")
    image = job.read_image()

    job.report(0.05, "Loading model")
    import torch
    torch.set_num_threads(int(os.environ.get("VHD_SEG_THREADS", torch.get_num_threads())))
    from cellpose import models
    # Falls back to the CPU when no GPU is visible
    model = models.Cellpose(gpu=True, model_type=p["model"])

    diameter = job.um_to_px(p["diameter_um"])
    # Enough context for a cell cut by the tile edge to be segmented whole
    overlap = int(round(2 * diameter))
    min_area = 0.05 * diameter ** 2

    total = count_blocks(image.shape, p["block_size"])
    polygons = []
    for i, (core, (y0, y1, x0, x1)) in enumerate(iter_blocks(image.shape, p["block_size"], overlap)):
        job.report(0.05 + 0.9 * i / total, f"Segmenting tile {i + 1}/{total}")
        tile = image[y0:y1, x0:x1]
        if tile.std() < 1:  # blank background
            continue
        masks = model.eval(
            tile, diameter=diameter, channels=[0, 0],
            flow_threshold=p["flow_threshold"], cellprob_threshold=p["cellprob_threshold"],
        )[0]
        polygons += labels_to_polygons(masks, offset=(y0, x0), core=core, min_area=min_area)

    job.report(0.97, "Writing polygons")
    job.write_polygons(polygons)


if __name__ == "__main__":
    main()
