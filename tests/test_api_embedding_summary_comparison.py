"""Embeddings, the dataset summary and sample comparison through the HTTP API."""
import json
import time
import unittest
from types import SimpleNamespace

from fastapi.testclient import TestClient

from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from tests.test_clustering_markers import make_layer
from tests.test_sample_comparison import _samples

WIDGETS = r'"name": "(Div|Paragraph|Button)"'      # the page loads no BokehJS widgets bundle


class APITestCase(unittest.TestCase):
    def open(self, samples):
        self.registry = SessionRegistry(ttl=3600)
        self.client = TestClient(create_app(self.registry, sweep_seconds=3600))
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)
        self.client.get("/api/load")
        self.session = self.registry.get(self.client.cookies[COOKIE])
        self.session.loaded = SimpleNamespace(samples=samples)
        self.session.views = build_views(self.session.loaded)

    def wait(self, path, key="job"):
        for _ in range(1200):
            job = self.client.get(path).json()[key]
            if job["state"] != "running":
                return job
            time.sleep(.05)
        self.fail(f"{path} still running")


class EmbeddingAndSummaryTests(APITestCase):
    def setUp(self):
        self.sample = make_layer().sample
        self.open({"bins": ("bins", self.sample, "Bins")})
        self.client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=181, y0=-1, y1=2))

    def test_embedding_figure_colours_hover_and_selection(self):
        self.assertEqual(self.client.get("/api/layers/bins/embedding/figure").status_code, 404)
        self.assertEqual(self.client.post("/api/layers/bins/embedding", json=dict(method="nope")).status_code, 422)
        self.assertEqual(self.client.post("/api/layers/bins/embedding", json=dict(method="pca")).status_code, 202)
        job = self.wait("/api/layers/bins/embedding")
        self.assertEqual(job, dict(name="embedding", state="done", message="Loaded 180 bin embeddings."))

        item = self.client.get("/api/layers/bins/embedding/figure", params=dict(colour_by="g0")).json()
        self.assertIn("doc", item)
        self.assertNotRegex(json.dumps(item), WIDGETS)
        recoloured = self.client.get("/api/layers/bins/embedding/colours", params=dict(colour_by="g0")).json()
        self.assertEqual(len(recoloured["codes"]), 180)
        self.assertGreater(len(recoloured["palette"]), 10)                   # a viridis ramp
        self.assertGreater(recoloured["version"], item["version"])

        names = self.client.post("/api/layers/bins/embedding/names",
                                 json=dict(row_version=recoloured["version"], rows=[0, 5])).json()["rows"]
        self.assertTrue(names[0][1].startswith("b0  ·  "))
        stale = self.client.post("/api/layers/bins/embedding/names", json=dict(row_version=0, rows=[0])).json()
        self.assertEqual(stale["rows"], [])

        selected = self.client.post("/api/layers/bins/embedding/selection", json=dict(indices=[0, 1, 2])).json()
        self.assertEqual(selected["selected"], 3)
        self.assertEqual(self.client.get("/api/layers/bins/embedding/selection").json()["indices"], [0, 1, 2])

        # Categories keep the clusters' assigned colours
        self.client.post("/api/layers/bins/clustering", json=dict(algorithm="leiden"))
        self.wait("/api/layers/bins/clusters")
        column = self.client.get("/api/layers/bins/embedding").json()["default_colour_by"]
        self.assertTrue(column)
        palette = self.client.get("/api/layers/bins/embedding/colours", params=dict(colour_by=column)).json()["palette"]
        clusters = self.client.get("/api/layers/bins/clusters").json()["clusters"]
        self.assertEqual(set(palette), {cluster["colour"] for cluster in clusters})

    def test_summary_qc_and_markers(self):
        stats, = self.client.get("/api/summary").json()["layers"]
        self.assertEqual((stats["id"], stats["observations"]), ("bins", 180))
        self.assertNotIn("not_sequenced", stats)
        self.assertAlmostEqual(stats["mean_umis"], round(float(self.sample.adata.X.sum(axis=1).mean()), 1))
        item = self.client.get("/api/layers/bins/qc").json()
        self.assertEqual((item["shown"], item["total"]), (180, 180))
        self.assertNotRegex(json.dumps(item), WIDGETS)
        self.assertEqual(self.client.get("/api/layers/bins/markers").json(), dict(clustered=False, markers=[]))
        self.client.post("/api/layers/bins/clustering", json=dict(algorithm="leiden"))
        self.wait("/api/layers/bins/clusters")
        markers = self.client.get("/api/layers/bins/markers").json()["markers"]
        self.assertTrue(markers and {"cluster", "gene", "score", "pval", "logfc"} <= set(markers[0]))


class ComparisonTests(APITestCase):
    def setUp(self):
        self.spots, self.bins = _samples()
        # The bins are given first: the spots (larger) are still the low-res reference
        self.open({"bins": ("spots", self.bins, "Bins"), "spots": ("spots", self.spots, "Spots")})

    def compare(self, **body):
        started = self.client.post("/api/comparison", json=dict(dict(sample1="bins", sample2="spots", min_coverage=.5),
                                                                **body))
        self.assertEqual(started.status_code, 202, started.text)
        return self.wait("/api/comparison")

    def test_options_run_results_plots_and_selection(self):
        options = self.client.get("/api/comparison").json()
        self.assertEqual((options["default_sample1"], options["default_sample2"]), ("spots", "bins"))
        spots = next(layer for layer in options["layers"] if layer["id"] == "spots")
        self.assertIn("spot_cluster", spots["clusters"])
        self.assertIn(dict(value="spot_deconv", label="spot_deconv (proportions)"), spots["celltypes"])
        self.assertEqual(self.client.post("/api/comparison", json=dict(sample1="bins", sample2="bins")).status_code, 422)
        self.assertEqual(self.client.get("/api/comparison/result").status_code, 404)

        job = self.compare(predicate="contains", sample1_cluster_key="bin_cluster", sample2_cluster_key="spot_cluster")
        self.assertEqual(job["state"], "done", job)
        self.assertIn("Compared 9 spots with 36 spots", job["message"])         # 2 x 2 bins inside each spot
        result = self.client.get("/api/comparison/result").json()
        self.assertEqual((result["low"], result["high"]), ("spots", "bins"))
        self.assertEqual(result["n_local"], 9)
        metrics = {row["metric"]: row for row in result["global_metrics"]}
        self.assertIn("cluster_ari", metrics)
        self.assertTrue(metrics["pseudobulk_pearson"]["help"])
        self.assertIn("pearson_hvg", result["metrics"])
        self.assertEqual(result["default_metric"], "pearson_hvg")
        self.assertEqual(len(result["genes"]), 12)
        self.assertEqual(len(result["genes"][0]), 1 + len(result["gene_columns"]))
        # Local metrics are columns of the low-res layer
        self.assertIn("cmp_pearson_hvg", self.spots.adata.obs)
        self.assertNotIn("cmp_pearson_hvg", self.bins.adata.obs)
        self.assertEqual(self.client.patch("/api/layers/spots/style",
                                           json=dict(colour_by="cmp_pearson_hvg")).status_code, 200)

        plots = self.client.get("/api/comparison/plots", params=dict(metric="pearson_hvg", x=result["default_x"],
                                                                    y=result["default_y"])).json()
        self.assertIn("doc", plots)
        self.assertNotRegex(json.dumps(plots), WIDGETS)
        self.assertEqual(self.client.get("/api/comparison/plots", params=dict(metric="nope", x="nope", y="nope"))
                         .status_code, 422)
        picked = self.client.post("/api/comparison/selection",
                                  json=dict(x=result["default_x"], y=result["default_y"], indices=[0, 1])).json()
        self.assertEqual((picked["layer"], picked["selected"]), ("spots", 2))

        # A new run replaces the metric columns (and the layer showing one re-reads it)
        self.assertEqual(self.compare(min_coverage=.99)["state"], "done")
        self.assertIn("cmp_pearson_hvg", self.spots.adata.obs)
        self.assertNotIn("cmp_cluster_ari", self.spots.adata.obs)


if __name__ == "__main__":
    unittest.main()
