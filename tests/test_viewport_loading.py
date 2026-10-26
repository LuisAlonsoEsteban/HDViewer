import threading
import time
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from shapely.geometry import box

from bokeh_server.components.glyphs import Patches
from src.sample import Cells
from tests.test_cell_loading import n_loaded, zoom_to


class QueueDoc:
    """Next-tick callbacks wait until `run_ticks` (like the server loop would)."""

    def __init__(self):
        self.callbacks = []
        self.lock = threading.Lock()

    def add_next_tick_callback(self, callback):
        with self.lock:
            self.callbacks.append(callback)

    def run_ticks(self):
        with self.lock:
            callbacks, self.callbacks = self.callbacks, []
        for callback in callbacks:
            callback()
        return len(callbacks)


def settle(doc, manager, timeout=5.0):
    """Run ticks until no load is running or queued."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        ran = doc.run_ticks()
        if not ran and not manager._in_flight and not doc.callbacks:
            return
        time.sleep(0.005)
    raise AssertionError("background load did not finish")


class SlowCells(Cells):
    """Counts fetches and can hold them until released."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_fetches = 0
        self.release = threading.Event()
        self.release.set()

    def fetch_data(self, *args, **kwargs):
        self.n_fetches += 1
        self.release.wait(5)
        return super().fetch_data(*args, **kwargs)


def make_layer(background=True, n=20):
    names = [f"c{i}_{j}" for i in range(n) for j in range(n)]
    gdf = gpd.GeoDataFrame({"index": names},
                           geometry=[box(50 * i, 50 * j, 50 * i + 10, 50 * j + 10)
                                     for i in range(n) for j in range(n)], index=names)
    adata = ad.AnnData(np.ones((n * n, 2)), obs=pd.DataFrame(index=names),
                       var=pd.DataFrame(index=["g0", "g1"]))
    plot = figure(x_range=(0, 1000), y_range=(0, 1000))
    layer = Patches(plot, SlowCells(gdf, adata), lazy_loading=True, alias="Cell", key="cells_vp")
    layer.add_parameter(attr="alpha", colname="g0", colname_hovered=True)
    layer.add_parameter(attr="color", static=True, val="#6baed6", not_sequenced_val="#b8b8b8")
    doc = QueueDoc()
    layer.attach_loader(doc)
    layer.loader.add_range_callback(layer, doc)
    layer.loader.background = background
    settle(doc, layer.manager)                 # initial (synchronous) load
    return layer, plot, doc


class ViewportLoadingTests(unittest.TestCase):
    def test_one_update_per_pan(self):
        layer, plot, doc = make_layer(background=False)
        layer.sample.n_fetches = 0
        zoom_to(plot, 100, 200)                # 4 range properties change
        self.assertEqual(len(doc.callbacks), 1)
        doc.run_ticks()
        self.assertEqual(layer.sample.n_fetches, 1)

    def test_background_load_adds_rows_on_a_later_tick(self):
        layer, plot, doc = make_layer()
        zoom_to(plot, 0, 90)
        doc.run_ticks()                        # update_bbox: clears and queues the load
        self.assertEqual(n_loaded(layer), 0)   # not added in the same tick
        settle(doc, layer.manager)
        self.assertGreater(n_loaded(layer), 0)
        self.assertLess(n_loaded(layer), 100)

    def test_stale_result_is_dropped(self):
        layer, plot, doc = make_layer()
        layer.sample.release.clear()           # hold the worker
        zoom_to(plot, 0, 90)
        doc.run_ticks()
        layer.manager.clear_sources()          # e.g. a new level or colour meanwhile
        layer.sample.release.set()
        settle(doc, layer.manager)
        self.assertEqual(n_loaded(layer), 0)

    def test_requests_during_a_load_are_merged(self):
        layer, plot, doc = make_layer()
        layer.sample.n_fetches = 0
        layer.sample.release.clear()
        zoom_to(plot, 0, 90)
        doc.run_ticks()
        for _ in range(5):                     # five more loads asked while one runs
            layer.manager.fetch_data_async(doc)
        layer.sample.release.set()
        settle(doc, layer.manager)
        self.assertEqual(layer.sample.n_fetches, 2)
        self.assertGreater(n_loaded(layer), 0)


if __name__ == "__main__":
    unittest.main()
