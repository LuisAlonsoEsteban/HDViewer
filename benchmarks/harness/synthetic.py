"""Synthetic Visium HD datasets in Space Ranger's output layout.

Every tool reads the same files through its usual reader, so the synthetic
sweep measures how each tool scales with the number of bins and cells:

    <out>/
      binned_outputs/square_{002,008,016}um/
        filtered_feature_bc_matrix.h5          (10x HDF5, in-tissue bins)
        spatial/tissue_positions.parquet       (all grid bins)
        spatial/scalefactors_json.json
        spatial/tissue_{hires,lowres}_image.png
      spatial/                                 (copy of the 2 µm spatial images)
      synthetic_tissue_image.tif               (full-resolution H&E-like image)
      cells.geojson                            (cell outlines, full-res pixels)
      meta.json

Expression has spatial domains with domain-specific gene profiles, so
clustering and differential expression find real structure. The 8 and 16 µm
levels are exact sums of the 2 µm counts, as in Space Ranger.

Usage: python -m benchmarks.harness.synthetic --bins 100000 --out data/benchmark_synthetic/100k
"""
import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

MICRONS_PER_PIXEL = 0.2738242950835738      # as the Mouse Brain slide
BIN_UM = 2
TISSUE_FRACTION = 0.7                       # of the square grid that is tissue
CELL_SPACING_UM = 12.0
N_DOMAINS = 8


def _tissue_mask(rows, cols, side, rng):
    """An irregular blob: an ellipse with a wavy border."""
    cy = cx = (side - 1) / 2
    dy, dx = (rows - cy) / (side / 2), (cols - cx) / (side / 2)
    angle = np.arctan2(dy, dx)
    phases = rng.uniform(0, 2 * np.pi, 3)
    border = 0.93 + 0.04 * np.sin(3 * angle + phases[0]) + 0.03 * np.sin(5 * angle + phases[1])
    return (dx / 1.0) ** 2 + (dy / 0.97) ** 2 <= border ** 2


def _domains(rows, cols, side, rng):
    """Domain label per bin: warped Voronoi regions of random seeds."""
    seeds = rng.uniform(0.1, 0.9, (N_DOMAINS, 2)) * side
    warp = side / 12
    wr = rows + warp * np.sin(cols / side * 2 * np.pi * 1.5)
    wc = cols + warp * np.sin(rows / side * 2 * np.pi * 2.0)
    labels = np.zeros(len(rows), dtype=np.int16)
    best = np.full(len(rows), np.inf)
    for k, (sr, sc_) in enumerate(seeds):
        d = (wr - sr) ** 2 + (wc - sc_) ** 2
        better = d < best
        labels[better], best[better] = k, d[better]
    return labels


def _counts(labels, n_genes, mean_umis, rng):
    """CSR (bins x genes): multinomial draws from each domain's gene profile."""
    base = rng.gamma(0.6, 1.0, n_genes)
    profiles = []
    for k in range(N_DOMAINS):
        p = base.copy()
        markers = rng.choice(n_genes, max(5, n_genes // 25), replace=False)
        p[markers] *= rng.uniform(4, 12, len(markers))
        p[k] *= 20                                   # Gene000k marks domain k
        profiles.append(p / p.sum())
    n = len(labels)
    totals = rng.poisson(rng.lognormal(np.log(mean_umis), 0.5, n))
    rows_out, cols_out = [], []
    for k in range(N_DOMAINS):
        idx = np.flatnonzero(labels == k)
        t = totals[idx]
        genes = rng.choice(n_genes, int(t.sum()), p=profiles[k]).astype(np.int32)
        rows_out.append(np.repeat(idx.astype(np.int64), t))
        cols_out.append(genes)
    r, c = np.concatenate(rows_out), np.concatenate(cols_out)
    m = sp.csr_matrix((np.ones(len(r), dtype=np.int32), (r, c)), shape=(n, n_genes))
    m.sum_duplicates()
    return m


def _write_10x_h5(path, matrix_csr, barcodes, genes):
    """10x Genomics feature-barcode HDF5 (genes x barcodes CSC == bins x genes CSR)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"), compression="gzip")
        g.create_dataset("data", data=matrix_csr.data.astype(np.int32), compression="gzip")
        g.create_dataset("indices", data=matrix_csr.indices.astype(np.int64), compression="gzip")
        g.create_dataset("indptr", data=matrix_csr.indptr.astype(np.int64), compression="gzip")
        g.create_dataset("shape", data=np.array([matrix_csr.shape[1], matrix_csr.shape[0]], dtype=np.int32))
        feat = g.create_group("features")
        n = len(genes)
        feat.create_dataset("id", data=np.array([f"SYNTH{i:05d}" for i in range(n)], dtype="S"))
        feat.create_dataset("name", data=np.array(genes, dtype="S"))
        feat.create_dataset("feature_type", data=np.array(["Gene Expression"] * n, dtype="S"))
        feat.create_dataset("genome", data=np.array(["synthetic"] * n, dtype="S"))
        feat.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))
        f.attrs["library_ids"] = np.array([b"synthetic"])
        f.attrs["chemistry_description"] = "Visium HD v1"
        f.attrs["filetype"] = "matrix"
        f.attrs["version"] = 2


def _barcodes(size_um, rows, cols):
    return [f"s_{size_um:03d}um_{r:05d}_{c:05d}-1" for r, c in zip(rows, cols)]


def _write_level(out, size_um, rows, cols, in_tissue, px_row, px_col, counts, genes, hires, lowres, scalef):
    folder = out / "binned_outputs" / f"square_{size_um:03d}um"
    spatial = folder / "spatial"
    spatial.mkdir(parents=True, exist_ok=True)
    barcodes = _barcodes(size_um, rows, cols)
    pd.DataFrame(dict(barcode=barcodes, in_tissue=in_tissue.astype(np.int32), array_row=rows.astype(np.int32),
                      array_col=cols.astype(np.int32), pxl_row_in_fullres=px_row,
                      pxl_col_in_fullres=px_col)).to_parquet(spatial / "tissue_positions.parquet", index=False)
    (spatial / "scalefactors_json.json").write_text(json.dumps(dict(
        spot_diameter_fullres=size_um / MICRONS_PER_PIXEL, bin_size_um=float(size_um),
        microns_per_pixel=MICRONS_PER_PIXEL, regist_target_img_scalef=scalef[0],
        tissue_hires_scalef=scalef[0], tissue_lowres_scalef=scalef[1],
        fiducial_diameter_fullres=size_um / MICRONS_PER_PIXEL * 50), indent=2))
    from PIL import Image
    Image.fromarray(hires).save(spatial / "tissue_hires_image.png")
    Image.fromarray(lowres).save(spatial / "tissue_lowres_image.png")
    keep = np.flatnonzero(in_tissue)
    _write_10x_h5(folder / "filtered_feature_bc_matrix.h5", counts, [barcodes[i] for i in keep], genes)


def _aggregate(rows, cols, in_tissue, counts_tissue, factor):
    """Sum 2 µm bins into factor x factor blocks (Space Ranger's 8/16 µm levels)."""
    side_r, side_c = rows.max() // factor + 1, cols.max() // factor + 1
    brow, bcol = rows // factor, cols // factor
    block = brow * side_c + bcol
    gr, gc = np.divmod(np.arange(side_r * side_c), side_c)
    tissue_idx = np.flatnonzero(in_tissue)
    agg = sp.csr_matrix((np.ones(len(tissue_idx), dtype=np.int32), (block[tissue_idx], np.arange(len(tissue_idx)))),
                        shape=(side_r * side_c, len(tissue_idx)))
    summed = (agg @ counts_tissue).tocsr()
    block_in_tissue = np.asarray(agg.sum(axis=1)).ravel() > 0
    return gr, gc, block_in_tissue, summed[block_in_tissue]


def _cells(side_um, mask_fn, rng):
    """Jittered-grid elliptical cells; returns centres (µm) and polygons (µm)."""
    n = int(side_um / CELL_SPACING_UM)
    g = (np.arange(n) + 0.5) * CELL_SPACING_UM
    cy, cx = [a.ravel() for a in np.meshgrid(g, g, indexing="ij")]
    cy = cy + rng.uniform(-2, 2, cy.size)
    cx = cx + rng.uniform(-2, 2, cx.size)
    keep = mask_fn(cy, cx)
    cy, cx = cy[keep], cx[keep]
    a = rng.uniform(3.5, 5.0, cy.size)
    b = a * rng.uniform(0.7, 1.0, cy.size)
    rot = rng.uniform(0, np.pi, cy.size)
    t = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    ex = a[:, None] * np.cos(t)[None]
    ey = b[:, None] * np.sin(t)[None]
    px = cx[:, None] + ex * np.cos(rot)[:, None] - ey * np.sin(rot)[:, None]
    py = cy[:, None] + ex * np.sin(rot)[:, None] + ey * np.cos(rot)[:, None]
    return cx, cy, px, py


def _write_image(path, size_px, origin_px, cells_px, rng, tile=2048):
    """H&E-like RGB TIFF (tiled), nuclei drawn at cell centres."""
    import cv2
    import tifffile

    cx, cy = cells_px
    order = np.argsort(cy)
    cx, cy = cx[order], cy[order]
    nucleus = max(2, int(round(2.5 / MICRONS_PER_PIXEL)))

    def tiles():
        for y0 in range(0, size_px, tile):
            for x0 in range(0, size_px, tile):
                h, w = min(tile, size_px - y0), min(tile, size_px - x0)
                img = np.empty((h, w, 3), dtype=np.uint8)
                noise = rng.normal(0, 6, (h, w)).astype(np.float32)
                img[..., 0] = np.clip(232 + noise, 0, 255)
                img[..., 1] = np.clip(196 + noise, 0, 255)
                img[..., 2] = np.clip(216 + noise, 0, 255)
                lo, hi = np.searchsorted(cy, [y0 - nucleus, y0 + h + nucleus])
                sel = (cx[lo:hi] >= x0 - nucleus) & (cx[lo:hi] < x0 + w + nucleus)
                for x, y in zip(cx[lo:hi][sel], cy[lo:hi][sel]):
                    cv2.circle(img, (int(x - x0), int(y - y0)), nucleus, (95, 60, 140), -1)
                yield img

    tifffile.imwrite(path, tiles(), shape=(size_px, size_px, 3), dtype=np.uint8, tile=(tile, tile),
                     photometric="rgb", compression="zlib", bigtiff=size_px * size_px * 3 > 3.5e9)


def generate(n_bins, out, n_genes=500, mean_umis=12.0, seed=0, with_image=True):
    rng = np.random.default_rng(seed)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    side = int(math.ceil(math.sqrt(n_bins / (TISSUE_FRACTION))))
    side = side + (-side % 8)                        # whole 16 µm blocks
    pitch_px = BIN_UM / MICRONS_PER_PIXEL
    margin_px = 200.0
    rows, cols = [a.ravel().astype(np.int64) for a in np.indices((side, side))]
    in_tissue = _tissue_mask(rows, cols, side, rng)
    px_row = margin_px + (rows + 0.5) * pitch_px
    px_col = margin_px + (cols + 0.5) * pitch_px
    labels = _domains(rows[in_tissue], cols[in_tissue], side, rng)
    genes = [f"Gene{i:04d}" for i in range(n_genes)]
    counts = _counts(labels, n_genes, mean_umis, rng)

    size_px = int(math.ceil(2 * margin_px + side * pitch_px))
    hires_scalef = min(1.0, 6000 / size_px)
    lowres_scalef = min(1.0, 600 / size_px)
    side_um = side * BIN_UM

    def tissue_at_um(y_um, x_um):
        r = np.clip((y_um / BIN_UM).astype(np.int64), 0, side - 1)
        c = np.clip((x_um / BIN_UM).astype(np.int64), 0, side - 1)
        return in_tissue[r * side + c]

    cx, cy, px, py = _cells(side_um, tissue_at_um, rng)
    to_px = lambda um: margin_px + um / MICRONS_PER_PIXEL
    cells_px = (to_px(cx), to_px(cy))

    image_path = out / "synthetic_tissue_image.tif"
    if with_image:
        _write_image(image_path, size_px, margin_px, cells_px, rng)
    import cv2
    full_small = _thumbnail(size_px, margin_px, pitch_px, side, in_tissue)
    hires = cv2.resize(full_small, (max(1, int(size_px * hires_scalef)),) * 2, interpolation=cv2.INTER_AREA)
    lowres = cv2.resize(full_small, (max(1, int(size_px * lowres_scalef)),) * 2, interpolation=cv2.INTER_AREA)

    levels = {}
    for size_um in (2, 8, 16):
        factor = size_um // BIN_UM
        if factor == 1:
            r, c, t, m = rows, cols, in_tissue, counts
        else:
            r, c, t, m = _aggregate(rows, cols, in_tissue, counts, factor)
        lp = size_um / MICRONS_PER_PIXEL
        _write_level(out, size_um, r, c, t, margin_px + (r + 0.5) * lp, margin_px + (c + 0.5) * lp,
                     m, genes, hires, lowres, (hires_scalef, lowres_scalef))
        levels[size_um] = dict(n_bins=int(t.sum()), n_grid=int(len(t)))
    # Space Ranger also writes the slide-level spatial folder
    spatial = out / "spatial"
    spatial.mkdir(exist_ok=True)
    from PIL import Image
    Image.fromarray(hires).save(spatial / "tissue_hires_image.png")
    Image.fromarray(lowres).save(spatial / "tissue_lowres_image.png")

    features = [dict(type="Feature", properties=dict(id=i + 1),
                     geometry=dict(type="Polygon", coordinates=[
                         [[round(float(x), 2), round(float(y), 2)] for x, y in
                          zip(to_px(px[i]), to_px(py[i]))] + [[round(float(to_px(px[i, 0])), 2),
                                                               round(float(to_px(py[i, 0])), 2)]]]))
                for i in range(len(cx))]
    with open(out / "cells.geojson", "w") as f:
        json.dump(dict(type="FeatureCollection", features=features), f)

    meta = dict(n_bins_2um=int(in_tissue.sum()), grid_side=side, image_px=size_px, n_genes=n_genes,
                n_cells=int(len(cx)), seed=seed, levels=levels, microns_per_pixel=MICRONS_PER_PIXEL,
                image=image_path.name if with_image else None)
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    return meta


def _thumbnail(size_px, margin_px, pitch_px, side, in_tissue, max_side=6000):
    """Low-resolution tissue picture for the hires/lowres PNGs."""
    import cv2
    scale = min(1.0, max_side / size_px)
    n = max(1, int(size_px * scale))
    img = np.full((n, n, 3), 245, dtype=np.uint8)
    mask = in_tissue.reshape(side, side).astype(np.uint8) * 255
    x0 = int(margin_px * scale)
    w = max(1, int(side * pitch_px * scale))
    mask = cv2.resize(mask, (w, w), interpolation=cv2.INTER_NEAREST)
    region = img[x0:x0 + w, x0:x0 + w]
    region[mask[:region.shape[0], :region.shape[1]] > 0] = (225, 180, 210)
    return img


def parse_size(text):
    text = str(text).lower()
    mult = {"k": 1_000, "m": 1_000_000}.get(text[-1], 1)
    return int(float(text[:-1] if text[-1] in "km" else text) * mult)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bins", required=True, help="in-tissue 2 µm bins, e.g. 100k, 2M")
    parser.add_argument("--out", required=True)
    parser.add_argument("--genes", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-image", action="store_true")
    args = parser.parse_args()
    print(json.dumps(generate(parse_size(args.bins), args.out, n_genes=args.genes, seed=args.seed,
                              with_image=not args.no_image), indent=2))
