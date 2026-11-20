"""Benchmark datasets: one spec per dataset, shared by every tool.

A spec holds the Space Ranger folder, the full-resolution image, the cell
outlines, the genes to colour by, the view boxes (harness/viewports.py) and
the analysis parameters. `prepare()` also builds a *reader mirror*: a folder of
symlinks laid out as spatialdata-io's Visium HD reader expects
(``<id>_tissue_image.<ext>``, ``<id>_feature_slice.h5``), so no file is written
into the original data folders.

The feature slice is a stub holding only the metadata the reader requires
(``hd_layout.file_format`` = "1.0"). Its transforms are only used to align the
CytAssist image, which the benchmark does not load.
"""
import json
import os
from pathlib import Path

import h5py

from . import viewports

REPO = Path(__file__).resolve().parents[2]
IDENTITY = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]


def _stub_feature_slice(path):
    metadata = dict(
        hd_layout_json=json.dumps(dict(file_format="1.0", transform=IDENTITY)),
        transform_matrices={k: IDENTITY for k in (
            "cytassist_colrow_to_spot_colrow", "spot_colrow_to_cytassist_colrow",
            "microscope_colrow_to_spot_colrow", "spot_colrow_to_microscope_colrow")},
        benchmark_stub=True)
    with h5py.File(path, "w") as f:
        f.attrs["metadata_json"] = json.dumps(metadata)


def reader_mirror(root, image_path, dataset_id, out):
    """Mirror `root` (files as symlinks) with spatialdata-io's expected names."""
    root, out = Path(root).resolve(), Path(out)
    for folder in ("binned_outputs", "spatial"):
        for source in (root / folder).rglob("*"):
            if source.is_dir() or source.suffix in (".gz", ".mtx", ".tsv"):
                continue
            target = out / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                os.symlink(source, target)
    image = Path(image_path).resolve()
    link = out / f"{dataset_id}_tissue_image{image.suffix}"
    if not link.exists():
        os.symlink(image, link)
    stub = out / f"{dataset_id}_feature_slice.h5"
    if not stub.exists():
        _stub_feature_slice(stub)
    return out


def microns_per_pixel(root, bin_size=2):
    path = Path(root) / "binned_outputs" / f"square_{bin_size:03d}um" / "spatial" / "scalefactors_json.json"
    scalefactors = json.loads(path.read_text())
    if "microns_per_pixel" in scalefactors:
        return float(scalefactors["microns_per_pixel"])
    return bin_size / float(scalefactors["spot_diameter_fullres"])


def prepare(name, entry, work_root, analysis):
    """The spec of dataset `name` (config entry), creating its reader mirror once."""
    root = (REPO / entry["root"]).resolve()
    image_path = (REPO / entry["image"]).resolve()
    cells_path = (REPO / entry["cells"]).resolve() if entry.get("cells") else None
    dataset_id = entry.get("dataset_id", name)
    mirror = reader_mirror(root, image_path, dataset_id, Path(work_root) / "mirrors" / name)
    mpp = microns_per_pixel(root)
    views = viewports.compute(root / "binned_outputs" / "square_002um" / "spatial" / "tissue_positions.parquet", mpp)
    return dict(dataset=name, dataset_id=dataset_id, root=str(root), sd_root=str(mirror),
                image_path=str(image_path), cells_path=str(cells_path) if cells_path else None,
                microns_per_pixel=mpp, genes=list(entry["genes"]), bin_sizes=list(entry["bin_sizes"]),
                cells_bin_size=2, viewports=views, analysis=dict(analysis), n_bins_2um=entry.get("n_bins_2um"))
