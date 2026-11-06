import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.startup.sample import Cells
from tests.fixtures import loaded, make_view, show


def make_cells(n=20):
    """n x n square cells of 10 units covering a 1000 x 1000 slide (every 50 units)."""
    names = [f"c{i}_{j}" for i in range(n) for j in range(n)]
    gdf = gpd.GeoDataFrame({"index": names},
                           geometry=[box(50 * i, 50 * j, 50 * i + 10, 50 * j + 10)
                                     for i in range(n) for j in range(n)], index=names)
    adata = ad.AnnData(np.ones((n * n, 2)), obs=pd.DataFrame(index=names),
                       var=pd.DataFrame(index=["g0", "g1"]))
    return Cells(gdf, adata)


def make_cell_view(n=20, sample=None):
    return make_view(make_cells(n) if sample is None else sample, kind="cells")


def zoom_to(view, x0, width, width_px=600, height_px=600):
    """Show the square [x0, x0 + width]² (as the browser's range change does)."""
    return show(view, x0, x0 + width, x0, x0 + width, width=width_px, height=height_px)


class CellLoadingTests(unittest.TestCase):
    def test_cells_load_at_every_zoom(self):
        view = make_cell_view()
        zoom_to(view, 0, 1000)
        self.assertEqual(len(loaded(view)), 400)              # whole slide: every cell

        zoom_to(view, 0, 90)
        self.assertGreater(len(loaded(view)), 0)
        self.assertLess(len(loaded(view)), 400 // 4)          # only the visible region

        zoom_to(view, 0, 1000)                                # back out: all again
        self.assertEqual(len(loaded(view)), 400)


if __name__ == "__main__":
    unittest.main()
