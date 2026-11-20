"""Shared view scenarios, in full-resolution image pixels (x = column, y = row).

Every tool works in (or converts to) full-resolution pixel coordinates: Space
Ranger's tissue positions, SpatialData's 'global' system and HDViewer's
SpatialTransform input. The orchestrator computes the boxes once per dataset
and passes them to every adapter, so all tools render the same regions.
"""
import numpy as np
import pandas as pd

# Square windows (µm) centred on the tissue, from overview to cell level
ZOOMS_UM = {"zoom_1000um": 1000.0, "zoom_250um": 250.0, "zoom_50um": 50.0}
PAN_STEPS = 10
PAN_FRACTION = 0.5          # each pan step moves half a window


def _box(cx, cy, half):
    return [float(cx - half), float(cx + half), float(cy - half), float(cy + half)]


def compute(positions_path, microns_per_pixel):
    """{'overview': box, 'zoom_*': box, 'pan': [box, ...]} with box = [x0, x1, y0, y1]."""
    pos = pd.read_parquet(positions_path, columns=["in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    pos = pos[pos["in_tissue"] == 1]
    x, y = pos["pxl_col_in_fullres"].to_numpy(), pos["pxl_row_in_fullres"].to_numpy()
    x0, x1, y0, y1 = np.percentile(x, 0.5), np.percentile(x, 99.5), np.percentile(y, 0.5), np.percentile(y, 99.5)
    side = max(x1 - x0, y1 - y0) * 1.05
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    # Centre zooms on tissue: the median bin position (never on background)
    mx, my = float(np.median(x)), float(np.median(y))
    views = {"overview": _box(cx, cy, side / 2)}
    for name, um in ZOOMS_UM.items():
        views[name] = _box(mx, my, um / microns_per_pixel / 2)
    # Pan around a 250 µm window: right, down, left, up (a closed loop over tissue)
    half = ZOOMS_UM["zoom_250um"] / microns_per_pixel / 2
    step = 2 * half * PAN_FRACTION
    moves = [(1, 0), (1, 0), (0, 1), (0, 1), (-1, 0), (-1, 0), (-1, 0), (0, -1), (0, -1), (1, 0)][:PAN_STEPS]
    px, py, pan = mx, my, []
    for dx, dy in moves:
        px, py = px + dx * step, py + dy * step
        pan.append(_box(px, py, half))
    views["pan"] = pan
    return views
