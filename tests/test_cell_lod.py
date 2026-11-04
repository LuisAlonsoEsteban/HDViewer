import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from hd_viewer.startup.sample import Cells
from tests.test_cell_loading import ImmediateDoc, n_loaded, zoom_to
from bokeh.plotting import figure
from hd_viewer.board.components.layer import Patches


def make_round_cells(n=4, radius=20.0):
    """n x n round cells (64-vertex polygons) every 50 units."""
    names = [f"c{i}_{j}" for i in range(n) for j in range(n)]
    geoms = [Point(50 * i + 25, 50 * j + 25).buffer(radius, quad_segs=16)
             for i in range(n) for j in range(n)]
    gdf = gpd.GeoDataFrame({"index": names}, geometry=geoms, index=names)
    adata = ad.AnnData(np.ones((n * n, 1)), obs=pd.DataFrame(index=names),
                       var=pd.DataFrame(index=["g0"]))
    return Cells(gdf, adata)


def n_vertices(data):
    return sum(len(xs) for xs in data["x"])


class CellLevelOfDetailTests(unittest.TestCase):
    def test_level_follows_pixel_size(self):
        cells = make_round_cells()
        self.assertIsNone(cells.get_lod_level(0.5))           # finer than base: original
        self.assertIsNone(cells.get_lod_level(1.1))           # within hysteresis of native
        self.assertEqual(cells.get_lod_level(1.2), "lod:0")
        self.assertEqual(cells.get_lod_level(5.0), "lod:2")
        self.assertEqual(cells.get_lod_level(1e9), f"lod:{Cells.LOD_MAX_LEVEL}")

    def test_hysteresis_keeps_level_near_boundary(self):
        cells = make_round_cells()
        self.assertEqual(cells.get_lod_level(4.2, current="lod:1"), "lod:1")   # just past 4
        self.assertEqual(cells.get_lod_level(1.9, current="lod:1"), "lod:1")
        self.assertEqual(cells.get_lod_level(6.0, current="lod:1"), "lod:2")
        self.assertIsNone(cells.get_lod_level(0.8, current="lod:0"))

    def test_simplified_coordinates(self):
        cells = make_round_cells()
        native = cells.fetch_data()
        coarse = cells.fetch_data(resolution_level="lod:2")
        self.assertEqual(list(native.index), list(coarse.index))
        self.assertEqual(list(native["obs_name"]), list(coarse["obs_name"]))
        self.assertEqual(len(native["x"].iloc[0]), 65)        # closed ring: 64 + 1
        self.assertLess(n_vertices(coarse), n_vertices(native) / 2)
        # Every cell stays a closed polygon, even at the coarsest level
        for xs, ys in zip(*cells.fetch_data(resolution_level="lod:8")[["x", "y"]].T.values):
            self.assertGreaterEqual(len(xs), 4)
            self.assertEqual((xs[0], ys[0]), (xs[-1], ys[-1]))

    def test_simplified_levels_are_computed_once_and_cropped(self):
        cells = make_round_cells()
        whole = cells.fetch_data(resolution_level="lod:2")
        self.assertEqual(list(cells._simplified), [4.0])
        cached = cells._simplified[4.0]

        part = cells.fetch_data(bbox=[0, 90, 0, 90], resolution_level="lod:2")
        self.assertIs(cells._simplified[4.0], cached)           # reused, not recomputed
        self.assertLess(len(part), len(whole))
        for name in part.index:                                 # same outline as uncropped
            self.assertEqual(part.loc[name, "x"], whole.loc[name, "x"])
            self.assertEqual(part.loc[name, "y"], whole.loc[name, "y"])

    def test_layer_reloads_with_a_coarser_level_when_zooming_out(self):
        plot = figure(x_range=(0, 200), y_range=(0, 200), width=100, height=100)
        layer = Patches(plot, make_round_cells(), lazy_loading=True, alias="Cell", key="cells_lod")
        layer.add_parameter(attr="alpha", colname="g0", colname_hovered=True)
        layer.add_parameter(attr="color", static=True, val="#6baed6", not_sequenced_val="#b8b8b8")
        layer.attach_loader(ImmediateDoc())
        layer.loader.add_range_callback(layer, ImmediateDoc())

        zoom_to(plot, 0, 50)                                  # 0.5 units per pixel
        self.assertIsNone(layer.manager.resolution_level)
        fine_per_cell = n_vertices(layer.manager.source.data) / n_loaded(layer)

        zoom_to(plot, 0, 200)                                 # 2 units per pixel
        self.assertEqual(layer.manager.resolution_level, "lod:1")
        self.assertEqual(n_loaded(layer), 16)
        self.assertLess(n_vertices(layer.manager.source.data) / 16, fine_per_cell / 2)


if __name__ == "__main__":
    unittest.main()
