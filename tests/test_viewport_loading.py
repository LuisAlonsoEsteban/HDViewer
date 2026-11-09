"""Viewport loading of LayerView: one fetch per needed update, none otherwise.

The Bokeh lazy loader's background loads (stale results, merged requests) are
gone: the API loads synchronously under the view's lock and rejects stale
viewport requests by sequence number (tests.test_api)."""
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.startup.sample import Cells
from tests.fixtures import loaded, make_view, show


class CountingCells(Cells):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_fetches = 0

    def fetch_data(self, *args, **kwargs):
        self.n_fetches += 1
        return super().fetch_data(*args, **kwargs)


def make_layer(n=20):
    names = [f"c{i}_{j}" for i in range(n) for j in range(n)]
    gdf = gpd.GeoDataFrame({"index": names},
                           geometry=[box(50 * i, 50 * j, 50 * i + 10, 50 * j + 10)
                                     for i in range(n) for j in range(n)], index=names)
    adata = ad.AnnData(np.ones((n * n, 2)), obs=pd.DataFrame(index=names),
                       var=pd.DataFrame(index=["g0", "g1"]))
    view = make_view(CountingCells(gdf, adata), kind="cells")
    show(view, 0, 1000, 0, 1000)
    return view


def pan_to(view, x0, width=200):
    return show(view, x0, x0 + width, x0, x0 + width, width=1000, height=1000)


class ViewportLoadingTests(unittest.TestCase):
    def test_one_fetch_per_pan(self):
        view = make_layer()
        view.sample.n_fetches = 0
        pan_to(view, 100)
        self.assertEqual(view.sample.n_fetches, 1)

    def test_small_pans_inside_the_static_box_fetch_nothing(self):
        view = make_layer()
        pan_to(view, 100)
        view.sample.n_fetches = 0
        before = len(loaded(view))
        self.assertEqual(pan_to(view, 110), [])          # centre still inside the static box
        self.assertEqual(view.sample.n_fetches, 0)
        self.assertEqual(len(loaded(view)), before)
        ops = pan_to(view, 400)                          # leaves it: one fetch, new rows appended
        self.assertEqual(view.sample.n_fetches, 1)
        self.assertIn("append", [op[0] for op in ops])
        self.assertGreater(len(loaded(view)), 0)

    def test_hidden_layers_fetch_nothing_and_reload_when_shown(self):
        view = make_layer()
        view.sample.n_fetches = 0
        view.set_visible(False)
        self.assertEqual(len(loaded(view)), 0)
        self.assertEqual(pan_to(view, 100), [])
        self.assertEqual(view.sample.n_fetches, 0)
        view.set_visible(True)
        pan_to(view, 100)
        self.assertEqual(view.sample.n_fetches, 1)
        self.assertGreater(len(loaded(view)), 0)

    def test_far_rows_are_pruned_once_over_budget(self):
        view = make_layer()
        pan_to(view, 0, width=100)                       # fine level: a corner only
        for x0 in range(100, 900, 100):                  # walk along the diagonal
            pan_to(view, x0, width=100)
            self.assertLessEqual(len(loaded(view)), 400)
        names = set(loaded(view))
        self.assertNotIn("c0_0", names)                  # the start corner was dropped
        self.assertIn("c17_17", names)
        self.assertEqual(len(view.rows.sink.data["x"]), len(loaded(view)))


if __name__ == "__main__":
    unittest.main()
