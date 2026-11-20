"""Helpers shared by the adapters (numpy only, so they import in every tool env)."""
import os


def level_dir(spec, bin_size):
    return os.path.join(spec["root"], "binned_outputs", f"square_{int(bin_size):03d}um")


def matrix_path(spec, bin_size=None):
    return os.path.join(level_dir(spec, bin_size or spec["bin_size"]), "filtered_feature_bc_matrix.h5")


def positions_path(spec, bin_size=None):
    return os.path.join(level_dir(spec, bin_size or spec["bin_size"]), "spatial", "tissue_positions.parquet")


def views(spec):
    """(name, box) of the single views: the overview, then the zooms (largest first)."""
    v = spec["viewports"]
    return [(name, v[name]) for name in v if name != "pan"]


def quadrant_groups(xy):
    """Tool-independent DE groups: the tissue's four quadrants around the median
    position (labels Q1..Q4), so every tool compares the same observations."""
    import numpy as np
    x, y = (np.asarray(a, dtype=float) for a in xy)
    mx, my = np.median(x), np.median(y)
    right, low = x >= mx, y >= my
    labels = np.where(right, np.where(low, "Q4", "Q1"), np.where(low, "Q3", "Q2"))
    return labels


def bin_area_px(spec, bin_size=None):
    """Side of a bin in full-resolution pixels."""
    return float(bin_size or spec["bin_size"]) / spec["microns_per_pixel"]
