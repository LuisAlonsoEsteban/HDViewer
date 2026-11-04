import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from shapely.geometry import box

from hd_viewer.board.components.layer import Patches
from hd_viewer.startup.sample import Cells


class ImmediateDoc:
    def add_next_tick_callback(self, callback):
        callback()


def make_cells(n=20):
    """n x n square cells of 10 units covering a 1000 x 1000 slide (every 50 units)."""
    names = [f"c{i}_{j}" for i in range(n) for j in range(n)]
    gdf = gpd.GeoDataFrame({"index": names},
                           geometry=[box(50 * i, 50 * j, 50 * i + 10, 50 * j + 10)
                                     for i in range(n) for j in range(n)], index=names)
    adata = ad.AnnData(np.ones((n * n, 2)), obs=pd.DataFrame(index=names),
                       var=pd.DataFrame(index=["g0", "g1"]))
    plot = figure(x_range=(0, 1000), y_range=(0, 1000))
    layer = Patches(plot, Cells(gdf, adata), lazy_loading=True, alias="Cell", key="cells_test")
    layer.add_parameter(attr="alpha", colname="g0", colname_hovered=True)
    layer.add_parameter(attr="color", static=True, val="#6baed6", not_sequenced_val="#b8b8b8")
    layer.attach_loader(ImmediateDoc())
    layer.loader.add_range_callback(layer, ImmediateDoc())
    return layer, plot


def n_loaded(layer):
    return len(layer.manager.source.data["obs_name"])


def zoom_to(plot, x0, width):
    plot.x_range.update(start=x0, end=x0 + width)
    plot.y_range.update(start=x0, end=x0 + width)


class CellLoadingTests(unittest.TestCase):
    def test_cells_load_at_every_zoom(self):
        layer, plot = make_cells()
        self.assertEqual(n_loaded(layer), 400)                # whole slide: every cell

        zoom_to(plot, 0, 90)
        self.assertGreater(n_loaded(layer), 0)
        self.assertLess(n_loaded(layer), 400 // 4)            # only the visible region

        zoom_to(plot, 0, 1000)                                # back out: all again
        self.assertEqual(n_loaded(layer), 400)


if __name__ == "__main__":
    unittest.main()
