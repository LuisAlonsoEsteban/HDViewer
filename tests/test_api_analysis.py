"""Clustering and groups through the HTTP API (on a synthetic 180-bin layer)."""
import json
import time
import unittest
from types import SimpleNamespace

import numpy as np
from fastapi.testclient import TestClient

from hd_viewer.api import frames
from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from tests.test_clustering_markers import make_layer


class AnalysisApiTests(unittest.TestCase):
    def setUp(self):
        self.registry = SessionRegistry(ttl=3600)
        self.client = TestClient(create_app(self.registry, sweep_seconds=3600))
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)
        self.client.get("/api/load")                       # creates the session
        session = self.registry.get(self.client.cookies[COOKIE])
        sample = make_layer().sample
        session.loaded = SimpleNamespace(samples={"bins": ("bins", sample, "Bins")})
        session.views = build_views(session.loaded)
        self.session, self.sample = session, sample
        n = sample.adata.n_obs
        self.view_body = dict(seq=1, x0=-1, x1=n + 1, y0=-1, y1=2, width=900, height=900)

    def viewport(self, seq):
        response = self.client.post("/api/layers/bins/viewport", json={**self.view_body, "seq": seq})
        return frames.decode(response.content)

    def wait_job(self):
        for _ in range(600):
            job = self.client.get("/api/layers/bins/clusters").json()["job"]
            if job["state"] != "running":
                return job
            time.sleep(.05)
        self.fail("clustering did not finish")

    def test_clustering_job_colours_the_layer_lists_clusters_and_hides_them(self):
        self.viewport(1)
        options = self.client.get("/api/clustering/options").json()
        self.assertIn("leiden", options["algorithms"])
        with self.client.websocket_connect("/api/events") as socket:
            json.loads(socket.receive_text())
            started = self.client.post("/api/layers/bins/clustering", json=dict(algorithm="leiden", resolution=1.0))
            self.assertEqual(started.status_code, 202)
            seen = []
            while not seen or not (seen[-1]["type"] == "job" and seen[-1]["state"] != "running"):
                seen.append(json.loads(socket.receive_text()))
        self.assertIn("layer_stale", [event["type"] for event in seen])
        self.assertEqual(seen[-1]["state"], "done")
        result = self.client.get("/api/layers/bins/clusters").json()
        clusters = result["clusters"]
        self.assertGreater(len(clusters), 1)
        self.assertTrue(result["shown"])
        self.assertTrue(all(c["markers"] for c in clusters if c["size"] >= 2))
        self.assertEqual(sum(c["size"] for c in clusters), self.sample.adata.n_obs)

        header, ops, _ = frames.decode(self.client.get("/api/layers/bins/state").content)
        column = result["column"]
        self.assertIn(column, ops[0][1])
        self.assertEqual(header["palette"][1:len(clusters) + 1], [c["colour"] for c in clusters])

        first = clusters[0]["label"]
        header, _, extra = frames.decode(self.client.patch(
            "/api/layers/bins/clusters", json=dict(visible=[first])).content)
        self.assertTrue(header["filtered"])
        labels = np.asarray(ops[0][1][column], dtype=object)
        np.testing.assert_array_equal(extra["visible"].astype(bool), labels == first)
        header, _, _ = frames.decode(self.client.patch(
            "/api/layers/bins/clusters", json=dict(colours={first: "#000000"})).content)
        self.assertIn("#000000", header["palette"])

    def test_groups_from_conditions(self):
        self.viewport(1)
        rejected = self.client.post("/api/layers/bins/groups/check", json=dict(column="g0", op=">", value="x"))
        self.assertEqual(rejected.status_code, 422)
        hint = self.client.post("/api/layers/bins/groups/check", json=dict(column="g0", op=">", value="5")).json()
        self.assertIn("max=", hint["hint"])
        created = self.client.post("/api/layers/bins/groups", json=dict(
            name="High g0", conditions=[dict(column="g0", op=">", value="5")])).json()
        group = created["groups"][0]
        self.assertEqual(group["name"], "High g0")
        expected = int((np.asarray(self.sample.adata.X[:, 0]).ravel() > 5).sum())
        self.assertEqual(group["size"], expected)
        self.assertTrue(created["shown"])
        duplicate = self.client.post("/api/layers/bins/groups", json=dict(
            name="High g0", conditions=[dict(column="g0", op=">", value="1")]))
        self.assertEqual(duplicate.status_code, 422)

        header, ops, _ = frames.decode(self.client.get("/api/layers/bins/state").content)
        self.assertIn("custom_groups", ops[0][1])
        header, ops, extra = frames.decode(self.client.patch(
            "/api/layers/bins/groups/High g0", json=dict(visible=False)).content)
        self.assertTrue(header["filtered"])
        self.assertFalse(extra["visible"].any())

        remaining = self.client.delete("/api/layers/bins/groups/High g0").json()
        self.assertEqual(remaining["groups"], [])
        self.assertFalse(remaining["shown"])


if __name__ == "__main__":
    unittest.main()
