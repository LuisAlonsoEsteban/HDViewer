"""HTTP API of the React front end, end to end on a tiny dataset."""
import json
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from fastapi.testclient import TestClient

from hd_viewer import settings
from hd_viewer.api import frames
from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.startup import datasets
from tests.test_core_pipeline import write_dataset


class ApiTestCase(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        root = Path(self.folder.name)
        (root / "data").mkdir()
        patcher = patch.object(settings, "UPLOAD_DIR", root / "uploads")
        patcher.start()
        self.addCleanup(patcher.stop)
        examples = patch.dict(datasets.AVAILABLE_DATASETS, {"Tiny": write_dataset(root / "data")}, clear=True)
        examples.start()
        self.addCleanup(examples.stop)
        self.registry = SessionRegistry(ttl=3600)
        self.client = TestClient(create_app(self.registry, sweep_seconds=3600))
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)

    def load_example(self):
        response = self.client.post("/api/load", json=dict(source="example", dataset_id="Tiny", pseudospots=False))
        self.assertEqual(response.status_code, 202, response.text)
        for _ in range(200):
            state = self.client.get("/api/load").json()
            if state["state"] != "running":
                break
            time.sleep(.05)
        self.assertEqual(state["state"], "ready", state)
        return state


class ConfigTests(ApiTestCase):
    def test_config_and_datasets(self):
        config = self.client.get("/api/config").json()
        self.assertEqual(config["branding"]["name"], "VisiumVisualizer")
        self.assertIn("log1p", [scale["value"] for scale in config["expression_scales"]])
        self.assertEqual(self.client.get("/api/datasets").json()[0]["id"], "Tiny")
        self.assertEqual(self.client.get("/static/bokehjs/bokeh.min.js").status_code, 200)


class UploadTests(ApiTestCase):
    def upload(self, field, name, data, chunk=4):
        start = self.client.post("/api/uploads", json=dict(field=field, filename=name, size=len(data)))
        if start.status_code != 200:
            return start
        upload_id = start.json()["upload_id"]
        for index in range(0, len(data), chunk):
            response = self.client.put(f"/api/uploads/{field}/{upload_id}/{index // chunk}",
                                       content=data[index:index + chunk])
        return response

    def test_chunks_are_assembled_in_the_session_folder_with_a_safe_name(self):
        response = self.upload("expression", "../../escape.h5", b"0123456789")
        self.assertTrue(response.json()["complete"])
        sid = self.client.cookies[COOKIE]
        session = self.registry.get(sid)
        path = session.uploads["expression"].path
        self.assertEqual(path, session.directory / "expression" / "escape.h5")
        self.assertEqual(path.read_bytes(), b"0123456789")
        self.assertEqual(self.client.get("/api/uploads").json()["expression"]["status"], "complete")

    def test_rejections(self):
        self.assertEqual(self.upload("expression", "counts.csv", b"x").status_code, 422)
        with patch.object(settings, "UPLOAD_MAX_BYTES", 4):
            self.assertEqual(self.upload("image", "image.png", b"0123456789").status_code, 422)
        start = self.client.post("/api/uploads", json=dict(field="positions", filename="p.parquet", size=8)).json()
        out_of_order = self.client.put(f"/api/uploads/positions/{start['upload_id']}/1", content=b"0123")
        self.assertEqual(out_of_order.status_code, 422)
        self.assertEqual(self.client.post("/api/load", json=dict(source="upload", bins_size=8)).status_code, 422)


class WorkspaceTests(ApiTestCase):
    def test_load_progress_figures_viewport_style_names_and_selection(self):
        self.assertEqual(self.client.get("/api/load").json(), {})        # creates the session
        with self.client.websocket_connect("/api/events") as socket:
            self.assertEqual(json.loads(socket.receive_text())["type"], "status")
            self.load_example()
            events = []
            while not events or events[-1]["type"] != "ready":
                events.append(json.loads(socket.receive_text()))
        self.assertIn(("tiles", 100), [(e.get("task"), e.get("pct")) for e in events])

        workspace = self.client.get("/api/workspace").json()
        layer = workspace["layers"][0]
        self.assertEqual((layer["id"], layer["kind"]), ("bins", "bins"))
        self.assertEqual(len(workspace["image"]["inverse"]), 9)

        item = self.client.get("/api/figures/map").json()
        names = {ref.get("attributes", {}).get("name") for ref in _objects(item)}
        for suffix in ("source", "mapper", "filter", "renderer", "hover"):
            self.assertIn(f"bins:{suffix}", names)
        self.assertIn("viewport_box", {ref.get("attributes", {}).get("name")
                                       for ref in _objects(self.client.get("/api/figures/thumbnail").json())})

        size = workspace["image"]["size"]
        response = self.client.post("/api/layers/bins/viewport",
                                    json=dict(seq=1, x0=0, x1=size, y0=0, y1=size, width=600, height=600))
        self.assertEqual(response.headers["content-type"], "application/octet-stream")
        header, ops, _ = frames.decode(response.content)
        self.assertEqual((header["base"], header["version"], header["rows"]), (0, 1, 36))  # 3 sequenced
        appended = [op for op in ops if op[0] == "append"][0][1]
        self.assertEqual(appended["x"].dtype, np.float32)
        self.assertIn("width", header["constants"])
        stale = self.client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=0, x1=size, y0=0, y1=size))
        self.assertEqual(stale.status_code, 409)

        gene = self.client.get("/api/layers/bins/columns", params=dict(q="")).json()[1]
        header, ops, _ = frames.decode(self.client.patch("/api/layers/bins/style",
                                                         json=dict(colour_by=gene, scale="log1p")).content)
        self.assertEqual(header["base"], 1)
        self.assertTrue(any(name == gene for op in ops for name in (op[2:3] and [op[1]] or op[1])))
        self.assertIn("(log-normalized)", header["tooltips"][-1][0])

        version = header["row_version"]
        named = dict(self.client.post("/api/layers/bins/names",
                                      json=dict(row_version=version, rows=list(range(36)))).json()["rows"])
        self.assertEqual(len(named), 36)
        self.assertEqual(self.client.post("/api/layers/bins/names", json=dict(row_version=-1, rows=[0])).json()["rows"], [])
        sequenced = next(row for row, name in named.items() if name in ("a", "b", "c"))
        self.assertEqual(self.client.post("/api/layers/bins/selection",
                                          json=dict(row_version=version, rows=[sequenced])).json()["selected"], 1)
        self.assertEqual(self.client.get("/api/layers/bins/selection").json()["rows"], [sequenced])

        full, ops, _ = frames.decode(self.client.get("/api/layers/bins/state").content)
        self.assertEqual(ops[0][0], "replace")
        self.assertEqual(len(ops[0][1]["x"]), 36)

    def test_tiles_are_cached_forever_and_misses_never(self):
        self.load_example()
        url = self.client.get("/api/workspace").json()["tiles"]["template"]
        url = url.replace("{Z}", "0").replace("{X}", "0").replace("{Y}", "0")
        hit = self.client.get(url)
        self.assertEqual(hit.status_code, 200)
        self.assertIn("immutable", hit.headers["cache-control"])
        self.assertEqual(self.client.get(url, headers={"if-none-match": hit.headers["etag"]}).status_code, 304)
        miss = self.client.get(url.replace("/0/0/0.", "/9/9/9."))
        self.assertEqual((miss.status_code, miss.headers["cache-control"]), (404, "no-store"))


class SessionTests(ApiTestCase):
    def test_sessions_are_private_and_closed_sessions_are_removed(self):
        self.load_example()
        other = TestClient(self.client.app)
        self.assertEqual(other.get("/api/workspace").status_code, 409)       # nothing loaded for them
        sid = self.client.cookies[COOKIE]
        session = self.registry.get(sid)
        session.directory.mkdir(parents=True, exist_ok=True)
        session.disconnected_at = time.monotonic() - 3600
        self.assertEqual(self.registry.sweep(), [sid])
        for _ in range(100):
            if not session.directory.exists():
                break
            time.sleep(.05)
        self.assertFalse(session.directory.exists())
        self.assertIsNone(self.registry.get(sid))


def _objects(value):
    if isinstance(value, dict):
        if value.get("type") == "object":
            yield value
        for item in value.values():
            yield from _objects(item)
    elif isinstance(value, list):
        for item in value:
            yield from _objects(item)


if __name__ == "__main__":
    unittest.main()
