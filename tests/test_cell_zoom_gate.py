import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from shapely.geometry import box

from bokeh_server.components.glyphs import Patches
from bokeh_server.components.zoom import ZoomBar
from src.sample import Cells

GATE = 100.0      # cells drawn when the view is at most 100 units wide (slide: 1000)


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
    layer.max_view_span = lambda: GATE
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


class CellZoomGateTests(unittest.TestCase):
    def test_cells_load_only_past_the_checkpoint(self):
        layer, plot = make_cells()
        self.assertEqual(n_loaded(layer), 0)                  # whole slide: nothing loaded
        self.assertIsNone(layer.manager.bbox)

        zoom_to(plot, 0, 90)                                  # narrower than the gate
        self.assertGreater(n_loaded(layer), 0)
        self.assertLess(n_loaded(layer), 400 // 4)            # only the visible region

        zoom_to(plot, 0, 500)                                 # back out: unloaded
        self.assertEqual(n_loaded(layer), 0)
        self.assertIsNone(layer.manager.bbox)

    def test_nothing_else_can_load_the_whole_slide_while_zoomed_out(self):
        layer, _ = make_cells()
        layer.manager.fetch_data()                            # no viewport: no-op
        layer.show_colname("g1")                              # recolour (used to fetch everything)
        layer.add_parameter(attr="line_width", static=True, val=2.0, not_sequenced_val=0.5)
        self.assertEqual(n_loaded(layer), 0)

    def test_zoom_bar(self):
        plot = figure(x_range=(0, 1000), y_range=(0, 1000))
        bar = ZoomBar(plot, full_span=1000)
        self.assertIn("Zoom <b style='color:#e5edf0'>1.0×</b>", bar.div.text)
        self.assertNotIn("cells", bar.div.text)               # no cells layer: no checkpoint

        bar.set_cells_gate(lambda: GATE, label="view ≤ 1 mm")
        self.assertEqual(bar.gate_zoom(), 10.0)
        self.assertIn("Cells hidden · zoom in to 10.0× (view ≤ 1 mm)", bar.div.text)
        self.assertFalse(bar.cells_visible())

        zoom_to(plot, 0, 50)                                  # 20x
        self.assertAlmostEqual(bar.zoom(), 20.0)
        self.assertIn("Cells shown", bar.div.text)
        self.assertTrue(bar.cells_visible())
        self.assertIn("left:47.", bar.div.text)               # checkpoint at log(10)/log(128)


if __name__ == "__main__":
    unittest.main()
