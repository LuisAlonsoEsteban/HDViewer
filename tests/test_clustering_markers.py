"""Marker genes of a clustering: a separate step, never stale, shown after the clusters."""
import threading
import time
import unittest
from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from fastapi.testclient import TestClient
from shapely.geometry import box

from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from hd_viewer.startup.sample import Bins
from tests.fixtures import make_view


def make_layer(n_per=60, n_types=3, n_genes=400):
    """The API view of a 180-bin layer (bin i = unit square [i, i+1] x [0, 1]):
    three cell types, each with its own block of 20 strongly expressed genes."""
    rng = np.random.default_rng(0)
    n = n_per * n_types
    X = rng.poisson(0.5, size=(n, n_genes)).astype(float)
    for t in range(n_types):
        X[t * n_per:(t + 1) * n_per, t * 20:(t + 1) * 20] += rng.poisson(15, size=(n_per, 20))
    names = [f"b{i}" for i in range(n)]
    index = pd.Index(names, name="barcode")
    gdf = gpd.GeoDataFrame({"index": names}, geometry=[box(i, 0, i + 1, 1) for i in range(n)], index=index)
    adata = ad.AnnData(X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    return make_view(Bins(gdf, adata, pixels_per_bin=1, bin_size=16), kind="bins", title="Bins")


def api_client(test, sample, kind="bins", layer_id="bins"):
    """(client, session) of a fresh API app whose session shows `sample` as one layer."""
    registry = SessionRegistry(ttl=3600)
    client = TestClient(create_app(registry, sweep_seconds=3600))
    client.__enter__()
    test.addCleanup(client.__exit__, None, None, None)
    client.get("/api/load")                       # creates the session
    session = registry.get(client.cookies[COOKIE])
    session.loaded = SimpleNamespace(samples={layer_id: (kind, sample, kind.title())})
    session.views = build_views(session.loaded)
    return client, session


class SampleMarkerTests(unittest.TestCase):
    def test_markers_are_a_separate_step_and_never_stale(self):
        clustering = make_layer().clustering
        clustering.run_clustering(resolution=1.0, compute_markers=False)
        self.assertEqual(clustering.get_marker_genes(), {})
        clustering.compute_markers()
        clusters = set(clustering.adata_processed.obs["leiden"].astype(str))
        self.assertEqual(set(clustering.get_marker_genes()), clusters)      # one entry per cluster

        clustering.run_clustering(resolution=1.5, compute_markers=False)   # new clustering
        self.assertEqual(clustering.get_marker_genes(), {})                # old ones hidden
        clustering.compute_markers()
        self.assertTrue(clustering.get_marker_genes())


class ClusteringJobTests(unittest.TestCase):
    def test_clusters_show_before_the_marker_genes(self):
        sample = make_layer().sample
        client, session = api_client(self, sample)
        client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=181, y0=-1, y1=2,
                                                           width=900, height=900))
        clustering = session.views["bins"].clustering
        release, ranking = threading.Event(), threading.Event()
        real_markers = clustering.compute_markers

        def slow_markers():
            ranking.set()
            release.wait(60)
            real_markers()
        clustering.compute_markers = slow_markers

        started = client.post("/api/layers/bins/clustering", json=dict(algorithm="leiden", resolution=1.0))
        self.assertEqual(started.status_code, 202)
        self.assertTrue(ranking.wait(120), client.get("/api/layers/bins/clusters").json()["job"])
        # Clustering done: clusters on the map and in the list, markers still running
        result = client.get("/api/layers/bins/clusters").json()
        self.assertTrue(clustering.is_clustered)
        self.assertTrue(result["shown"])
        self.assertEqual(session.views["bins"].colour_by["name"], clustering.cluster_col)
        self.assertTrue(result["clusters"])
        self.assertFalse(any(c["markers"] for c in result["clusters"]))
        self.assertEqual(result["job"]["state"], "running")
        self.assertIn("marker genes", result["job"]["message"])
        # Nothing else may run on the layer meanwhile
        self.assertEqual(client.post("/api/layers/bins/clustering", json=dict(resolution=1.0)).status_code, 409)

        release.set()
        deadline = time.time() + 120
        while client.get("/api/layers/bins/clusters").json()["job"]["state"] == "running" \
                and time.time() < deadline:
            time.sleep(0.2)
        result = client.get("/api/layers/bins/clusters").json()
        self.assertEqual(result["job"]["state"], "done")
        self.assertIn("marker genes ready", result["job"]["message"])
        self.assertTrue(any(gene.startswith("g") for c in result["clusters"] for gene in c["markers"]))


if __name__ == "__main__":
    unittest.main()
