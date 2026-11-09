"""Deconvolution and metadata comparison through the HTTP API (synthetic 180-bin layer).

Bin i is the unit square [i, i+1] x [0, 1]; bins 0-59, 60-119 and 120-179 are three
cell types, each with its own block of 20 strongly expressed genes; the reference has
one profile per type. NNLS runs for real with the local runner."""
import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from fastapi.testclient import TestClient

from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from hd_viewer.core import deconvolution
from tests.test_clustering_markers import make_layer


def reference(n_per=5, n_genes=400):
    rng = np.random.default_rng(1)
    X = rng.poisson(0.2, size=(3 * n_per, n_genes)).astype(float)
    for t in range(3):
        X[t * n_per:(t + 1) * n_per, t * 20:(t + 1) * 20] += 30
    obs = pd.DataFrame({"cell_type": np.repeat(["A", "B", "C"], n_per), "donor": "d1"},
                       index=[f"c{i}" for i in range(3 * n_per)])
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))


class DeconvolutionAndMetadataTests(unittest.TestCase):
    def setUp(self):
        environment = patch.dict(os.environ, {"VHD_DECONV_RUNNER": "local"})
        environment.start()
        self.addCleanup(environment.stop)
        deconvolution.default_service.cache_clear()
        self.addCleanup(deconvolution.default_service.cache_clear)
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.base = Path(folder.name)

        self.registry = SessionRegistry(ttl=3600)
        self.client = TestClient(create_app(self.registry, sweep_seconds=3600))
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)
        self.client.get("/api/load")
        self.session = self.registry.get(self.client.cookies[COOKIE])
        self.sample = make_layer().sample
        self.session.loaded = SimpleNamespace(samples={"bins": ("bins", self.sample, "Bins")}, reference=reference(),
                                              reference_name="ref.h5ad", reference_label_column="cell_type",
                                              base_dir=self.base, dataset=SimpleNamespace(bins_size=16))
        self.session.views = build_views(self.session.loaded)
        self.client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=181, y0=-1, y1=2))

    def deconvolve(self):
        started = self.client.post("/api/layers/bins/deconvolution", json=dict(method="nnls", label_column="cell_type"))
        self.assertEqual(started.status_code, 202, started.text)
        for _ in range(1200):
            state = self.client.get("/api/layers/bins/deconvolution").json()
            if state["job"]["state"] != "running":
                break
            time.sleep(.05)
        self.assertEqual(state["job"]["state"], "done", state["job"])
        return state

    def group(self, name, x0, x1):
        roi = self.client.post("/api/layers/bins/rois", json=dict(kind="box", name=name, box=[x0, .2, x1, .8])).json()
        self.client.post(f"/api/layers/bins/rois/{roi['rois'][-1]['id']}/group")
        return f"ROI: {name}"

    def test_methods_reference_and_validation(self):
        state = self.client.get("/api/layers/bins/deconvolution").json()
        self.assertEqual(state["unavailable"], "")
        self.assertEqual(state["reference"], dict(name="ref.h5ad", cells=15, genes=400))
        self.assertEqual(state["label_column"], "cell_type")
        self.assertEqual(state["columns"], ["cell_type", "donor"])
        nnls = next(method for method in state["methods"] if method["name"] == "nnls")
        self.assertEqual(nnls["unavailable"], "")
        self.assertEqual(nnls["params"][0]["name"], "batch_size")
        self.assertEqual(state["runs"], [])
        bad = self.client.post("/api/layers/bins/deconvolution/check", json=dict(method="nnls", label_column="nope"))
        self.assertEqual(bad.status_code, 422)
        self.assertIn("missing", bad.json()["detail"])
        good = self.client.post("/api/layers/bins/deconvolution/check", json=dict(method="nnls", label_column="cell_type"))
        self.assertEqual(good.status_code, 200)

        self.session.loaded.reference = None
        state = self.client.get("/api/layers/bins/deconvolution").json()
        self.assertIn("No single-cell reference", state["unavailable"])
        self.assertEqual(self.client.post("/api/layers/bins/deconvolution",
                                          json=dict(method="nnls", label_column="cell_type")).status_code, 422)

    def test_nnls_run_publishes_columns_cached_under_the_dataset(self):
        state = self.deconvolve()
        run, = state["runs"]
        self.assertEqual(run["cell_types"], ["A", "B", "C"])
        self.assertEqual(run["fraction_columns"], ["deconv_1:A", "deconv_1:B", "deconv_1:C"])
        obs = self.sample.adata.obs
        self.assertGreater(obs["deconv_1:A"].iloc[:60].mean(), .8)        # type-0 bins are mostly A
        self.assertGreater(obs["deconv_1:C"].iloc[120:].mean(), .8)
        self.assertIn("deconv_1:A", self.client.get("/api/layers/bins/columns", params=dict(q="deconv_1")).json())
        self.assertEqual(self.client.patch("/api/layers/bins/style", json=dict(colour_by="deconv_1:B")).status_code, 200)
        self.assertTrue(list((self.base / "deconvolution" / "square_016um" / "bins" / "nnls").iterdir()))
        self.assertEqual(self.deconvolve()["runs"][1]["id"], "2")             # cached result, new run

    def test_cancelling_kills_the_run_and_keeps_no_result(self):
        service = deconvolution.default_service()
        started, finished = [], []
        real_run = service.runner.run

        def slow(plugin, input_dir, output_dir, on_progress, log_path):
            started.append(log_path)
            try:
                with open(log_path, "w") as log:
                    service.runner._stream([sys.executable, "-c", "import time; print('go', flush=True); time.sleep(120)"],
                                           None, log, on_progress)
            finally:
                finished.append(time.monotonic())

        with patch.object(service.runner, "run", slow):
            self.assertEqual(self.client.post("/api/layers/bins/deconvolution/cancel").status_code, 409)
            self.assertEqual(self.client.post("/api/layers/bins/deconvolution",
                                              json=dict(method="nnls", label_column="cell_type")).status_code, 202)
            for _ in range(600):
                if started:
                    break
                time.sleep(.05)
            self.assertTrue(started, "the plugin never started")
            time.sleep(.5)                                                     # the process is running
            cancelled = self.client.post("/api/layers/bins/deconvolution/cancel")
            self.assertEqual(cancelled.status_code, 200, cancelled.text)
            self.assertEqual(cancelled.json()["state"], "cancelled")
            for _ in range(200):
                if finished:
                    break
                time.sleep(.05)
            self.assertTrue(finished, "the plugin process was not killed")
            time.sleep(.3)                                                     # the worker thread unwinds
            state = self.client.get("/api/layers/bins/deconvolution").json()
            self.assertEqual(state["job"]["state"], "cancelled")              # not overwritten as failed
            self.assertEqual(state["runs"], [])
            self.assertIn("cancelled by the user", Path(started[0]).read_text())
        self.assertIs(service.runner.run.__func__, real_run.__func__)
        self.assertEqual(self.deconvolve()["runs"][0]["id"], "1")              # a new run works afterwards

    def test_metadata_comparison_of_fractions_between_groups_and_against_the_rest(self):
        options = self.client.get("/api/layers/bins/metadata/options").json()
        self.assertEqual(options["default_preset"], "numeric")
        self.deconvolve()
        first, second = self.group("A", .5, 39.5), self.group("B", 60.5, 99.5)
        options = self.client.get("/api/layers/bins/metadata/options").json()
        self.assertEqual(options["default_preset"], "1:fraction")
        fractions = next(p for p in options["presets"] if p["key"] == "1:fraction")
        self.assertEqual([c["label"] for c in fractions["columns"]], ["A", "B", "C"])
        self.assertEqual([g["name"] for g in options["groups"]], [first, second])
        self.assertFalse(options["ready"])

        columns = [c["name"] for c in fractions["columns"]]
        one = self.client.post("/api/layers/bins/metadata/compare", json=dict(groups=[first], preset="1:fraction",
                                                                           columns=columns))
        self.assertEqual(one.status_code, 422)
        compared = self.client.post("/api/layers/bins/metadata/compare",
                                    json=dict(groups=[first, second], preset="1:fraction", columns=columns))
        self.assertEqual(compared.json(), dict(ready=True, plots=1))
        item = self.client.get("/api/layers/bins/metadata/plots").json()
        # The page loads no BokehJS widgets bundle: no Div & co. in served figures
        self.assertNotRegex(json.dumps(item), r'"name": "(Div|Paragraph|Button)"')
        rest = self.client.post("/api/layers/bins/metadata/compare",
                                json=dict(groups=[first, second], include_rest=True, preset="1:fraction",
                                          columns=columns, statistic="median"))
        self.assertEqual(rest.json()["plots"], 2)
        self.assertEqual(self.client.post("/api/layers/bins/metadata/compare",
                                          json=dict(groups=[first, second], preset="1:fraction",
                                                    columns=["g0"])).status_code, 422)

        # A new group changes The rest: the comparison is dropped
        self.group("C", 120.5, 150.5)
        self.assertFalse(self.client.get("/api/layers/bins/metadata/options").json()["ready"])
        self.client.post("/api/layers/bins/metadata/compare",
                         json=dict(groups=[first, second], preset="1:fraction", columns=columns))
        self.client.delete(f"/api/layers/bins/groups/{second}")
        self.assertFalse(self.client.get("/api/layers/bins/metadata/options").json()["ready"])
        self.assertEqual(self.client.get("/api/layers/bins/metadata/plots").status_code, 404)


if __name__ == "__main__":
    unittest.main()
