"""ROIs and differential expression through the HTTP API (synthetic 180-bin layer).

Bin i is the unit square [i, i+1] x [0, 1]; bins 0-59, 60-119 and 120-179 are three
cell types, each with its own block of 20 strongly expressed genes."""
import json
import time
import unittest
from types import SimpleNamespace

from fastapi.testclient import TestClient

from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from tests.test_clustering_markers import make_layer


class RegionsAndDETests(unittest.TestCase):
    def setUp(self):
        self.registry = SessionRegistry(ttl=3600)
        self.client = TestClient(create_app(self.registry, sweep_seconds=3600))
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)
        self.client.get("/api/load")
        session = self.registry.get(self.client.cookies[COOKIE])
        self.sample = make_layer().sample
        session.loaded = SimpleNamespace(samples={"bins": ("bins", self.sample, "Bins")})
        session.views = build_views(session.loaded)
        self.client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=181, y0=-1, y1=2))

    def roi(self, **shape):
        response = self.client.post("/api/layers/bins/rois", json=shape)
        self.assertEqual(response.status_code, 200, response.text)
        return response.json()["rois"][-1]

    def test_shapes_membership_groups_and_selection(self):
        box = self.roi(kind="box", name="First ten", box=[0.5, 0.2, 9.5, 0.8])
        self.assertEqual(box["count"], 10)
        self.assertTrue(box["xs"] and box["ys"])
        lasso = self.roi(kind="lasso", xs=[20.5, 24.5, 24.5], ys=[0.1, 0.1, 0.9])
        self.assertEqual(lasso["name"], "ROI 2")
        self.assertEqual(lasso["count"], 5)                  # bins 20-24 touch the triangle
        point = self.roi(kind="point", centre=[50.5, 0.5], radius=1.2)
        self.assertEqual(point["count"], 3)                  # bins 49-51 touch the circle
        self.assertEqual(self.client.post("/api/layers/bins/rois", json=dict(kind="box", box=[1, 1, 1, 1])).status_code, 422)
        self.assertEqual(self.client.post("/api/layers/bins/rois", json=dict(kind="box", name="First ten",
                                                                              box=[0, 0, 2, 1])).status_code, 422)

        self.assertIn("b100", self.client.get("/api/layers/bins/barcodes", params=dict(q="b10")).json())
        seeded = self.roi(kind="barcode", barcode="b100", radius=0.6)
        self.assertEqual(seeded["count"], 3)

        made = self.client.post(f"/api/layers/bins/rois/{box['id']}/group").json()
        self.assertEqual(made["group"], "ROI: First ten")
        self.assertEqual(made["groups"][0]["size"], 10)
        self.assertEqual(made["rois"][0]["group"], "ROI: First ten")
        selected = self.client.post(f"/api/layers/bins/rois/{box['id']}/select").json()
        self.assertEqual(selected["selected"], 10)
        self.assertEqual(len(selected["rows"]), 10)
        self.client.post(f"/api/layers/bins/rois/{lasso['id']}/select")
        remaining = self.client.post(f"/api/layers/bins/rois/{box['id']}/unselect").json()
        self.assertEqual(remaining["selected"], 5)            # the lasso's bins stay selected

        removed = self.client.delete(f"/api/layers/bins/rois/{box['id']}").json()
        self.assertEqual(len(removed["rois"]), 3)
        self.assertEqual(removed["groups"], [])               # its group went with it
        focus = self.client.post(f"/api/layers/bins/rois/{lasso['id']}/focus").json()
        self.assertLess(focus["x0"], 20.5)

    def test_group_seeded_roi_preview_and_save(self):
        source = self.roi(kind="box", name="Seeds", box=[10.2, 0.2, 12.8, 0.8])     # bins 10-12
        self.client.post(f"/api/layers/bins/rois/{source['id']}/group")
        shape = dict(kind="group", group="ROI: Seeds", radius=1.2, exclude_group=True)
        preview = self.client.post("/api/layers/bins/rois/preview", json=shape).json()
        self.assertEqual(preview["seeds"], 3)
        self.assertLess(preview["bounds"][0], 10)
        saved = self.roi(**shape)
        self.assertEqual(saved["seed_group"], "ROI: Seeds")
        self.assertTrue(saved["group_excluded"])
        self.assertEqual(saved["count"], 2)                   # bins 9 and 13: the reach is x 9.3-13.7
        missing = self.client.post("/api/layers/bins/rois/preview", json=dict(shape, group="Nope"))
        self.assertEqual(missing.status_code, 422)

    def run_de(self, **body):
        started = self.client.post("/api/layers/bins/de", json=body)
        self.assertEqual(started.status_code, 202, started.text)
        for _ in range(600):
            job = self.client.get("/api/layers/bins/clusters").json()["job"]
            if job["state"] != "running":
                break
            time.sleep(.05)
        self.assertEqual(job["state"], "done", job)

    def test_deleting_a_compared_group_drops_the_results(self):
        first = self.roi(kind="box", name="A", box=[0.5, 0.2, 39.5, 0.8])
        second = self.roi(kind="box", name="B", box=[60.5, 0.2, 99.5, 0.8])
        for roi in (first, second):
            self.client.post(f"/api/layers/bins/rois/{roi['id']}/group")
        self.run_de(groups=["ROI: A"])
        self.assertTrue(self.client.get("/api/layers/bins/de/options").json()["ready"])
        self.client.delete(f"/api/layers/bins/rois/{second['id']}")          # not compared: kept
        self.assertTrue(self.client.get("/api/layers/bins/de/options").json()["ready"])
        self.client.delete("/api/layers/bins/groups/ROI: A")
        self.assertFalse(self.client.get("/api/layers/bins/de/options").json()["ready"])
        self.assertEqual(self.client.get("/api/layers/bins/de/table").status_code, 404)

    def test_de_between_two_roi_groups(self):
        first = self.roi(kind="box", name="A", box=[0.5, 0.2, 39.5, 0.8])
        second = self.roi(kind="box", name="B", box=[60.5, 0.2, 99.5, 0.8])
        for roi in (first, second):
            self.client.post(f"/api/layers/bins/rois/{roi['id']}/group")
        options = self.client.get("/api/layers/bins/de/options").json()
        self.assertEqual([g["name"] for g in options["groups"]], ["ROI: A", "ROI: B"])
        self.assertEqual(self.client.post("/api/layers/bins/de", json=dict(groups=["ROI: A"], mode=1)).status_code, 422)
        self.run_de(groups=["ROI: A", "ROI: B"], mode=1)

        table = self.client.get("/api/layers/bins/de/table", params=dict(limit=5)).json()
        self.assertEqual(table["sets"], ["ROI: A", "ROI: B"])
        self.assertTrue(all(colour.startswith("#") for colour in table["colours"]))
        self.assertEqual(len(table["rows"]), 5)
        self.assertIn("ROI: A_lfc", table["columns"])
        top = table["rows"][0][0]
        self.assertTrue(top.startswith("g"))
        self.assertLess(int(top[1:]), 40)                     # a type-0 or type-1 marker gene
        found = self.client.get("/api/layers/bins/de/table", params=dict(q="g5")).json()
        self.assertTrue(all("g5" in row[0] for row in found["rows"]))
        ordered = self.client.get("/api/layers/bins/de/table", params=dict(sort="ROI: A_lfc", desc=True)).json()["rows"]
        values = [row[1] for row in ordered if row[1] is not None]
        self.assertEqual(values, sorted(values, reverse=True))
        plots = self.client.get("/api/layers/bins/de/plots").json()
        self.assertIn("doc", plots)
        self.assertNotRegex(json.dumps(plots), r'"name": "(Div|Paragraph|Button)"')   # no widgets bundle


if __name__ == "__main__":
    unittest.main()
