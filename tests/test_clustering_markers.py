import threading
import time
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.plotting import figure
from shapely.geometry import box

from bokeh_server.components.glyphs import Squares
from src.sample import Bins


def make_layer(n_per=60, n_types=3, n_genes=400):
    """Three cell types, each with its own block of 20 strongly expressed genes."""
    rng = np.random.default_rng(0)
    n = n_per * n_types
    X = rng.poisson(0.5, size=(n, n_genes)).astype(float)
    for t in range(n_types):
        X[t * n_per:(t + 1) * n_per, t * 20:(t + 1) * 20] += rng.poisson(15, size=(n_per, 20))
    names = [f"b{i}" for i in range(n)]
    index = pd.Index(names, name="barcode")
    gdf = gpd.GeoDataFrame({"index": names}, geometry=[box(i, 0, i + 1, 1) for i in range(n)], index=index)
    adata = ad.AnnData(X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    layer = Squares(figure(), Bins(gdf, adata, pixels_per_bin=1, bin_size=16), lazy_loading=True, alias="Bin")
    layer.add_parameter(attr="alpha", colname="g0", colname_hovered=True)
    layer.add_parameter(attr="color", static=True, val="#6baed6", not_sequenced_val="#b8b8b8")
    return layer


class ImmediateDoc(Document):
    """Runs next-tick callbacks at once (from whichever thread schedules them)."""
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


class SampleMarkerTests(unittest.TestCase):
    def test_markers_are_a_separate_step_and_never_stale(self):
        sample = make_layer().sample
        sample.run_clustering(resolution=0.5, compute_markers=False)
        self.assertEqual(sample.get_marker_genes(), {})
        sample.compute_marker_genes()
        clusters = set(sample.adata_processed.obs["leiden"].astype(str))
        self.assertEqual(set(sample.get_marker_genes()), clusters)      # one entry per cluster

        sample.run_clustering(resolution=1.5, compute_markers=False)   # new clustering
        self.assertEqual(sample.get_marker_genes(), {})                # old ones hidden
        sample.compute_marker_genes()
        self.assertTrue(sample.get_marker_genes())


class ClusteringTabTests(unittest.TestCase):
    def test_clusters_show_before_the_marker_genes(self):
        layer = make_layer()
        clustering = layer.clustering_manager
        release, ranking = threading.Event(), threading.Event()
        real_markers = clustering.compute_markers

        def slow_markers():
            ranking.set()
            release.wait(60)
            real_markers()
        clustering.compute_markers = slow_markers

        # As in the app: the Clustering tab lives in the layer's full widget panel
        tabs = layer.get_widgets(ImmediateDoc())
        from bokeh.models import Button, Div, Slider
        widgets = next(t for t in tabs.tabs if t.title == "Clustering").child.children
        slider = next(w for w in widgets if isinstance(w, Slider))
        button = next(w for w in widgets if isinstance(w, Button))
        status = widgets[widgets.index(button) + 1]            # the status line below it
        self.assertIsInstance(status, Div)
        slider.value = 0.5
        button._trigger_event(ButtonClick(button))

        self.assertTrue(ranking.wait(120), status.text)
        # Leiden done: clusters on the map and in the list, markers still running
        self.assertTrue(clustering.is_clustered)
        self.assertEqual(layer.manager.get_param("color").colname, clustering.cluster_col)
        self.assertIn("computing marker genes", str([w.text for w in clustering._cluster_layout.children
                                                     if hasattr(w, "text")]))
        self.assertIn("marker genes", status.text)
        self.assertTrue(button.disabled)

        release.set()
        deadline = time.time() + 120
        while "marker genes ready" not in status.text and time.time() < deadline:
            time.sleep(0.2)
        self.assertIn("marker genes ready", status.text)
        self.assertFalse(button.disabled)
        gene_lists = [w.text for w in clustering._cluster_layout.children if hasattr(w, "text")]
        self.assertNotIn("computing marker genes", str(gene_lists))
        self.assertTrue(any("g" in t for t in gene_lists))


if __name__ == "__main__":
    unittest.main()
