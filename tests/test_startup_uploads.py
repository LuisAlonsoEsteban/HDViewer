"""Uploads and the loading screen: chunk assembly (core.uploads.UploadField) and
when a load may start (POST /api/uploads, PUT chunks, POST /api/load)."""
import tempfile
import threading
import time
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from hd_viewer import settings
from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.core.session import session_upload_dir
from hd_viewer.core.uploads import UploadField
from hd_viewer.startup import datasets

REQUIRED = (("expression", "counts.h5"), ("positions", "positions.parquet"), ("image", "tissue.png"))


class UploadApiTestCase(unittest.TestCase):
    """An app whose only example is "Tiny" (as tests.test_api.ApiTestCase, which is not
    imported: tests.test_api imports tests.test_deconvolution, which imports this module).
    Drives the chunked-upload protocol like frontend/src (POST, then PUT chunks)."""

    def setUp(self):
        from tests.test_core_pipeline import write_dataset
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
        quiet = patch("hd_viewer.api.routers.load.traceback")  # failed loads are expected here
        quiet.start()
        self.addCleanup(quiet.stop)
        self.client.get("/api/uploads")                      # creates the session

    def load_example(self):
        response = self.post_load(source="example", dataset_id="Tiny")
        self.assertEqual(response.status_code, 202, response.text)
        state = self.wait_load()
        self.assertEqual(state["state"], "ready", state)
        return state

    @property
    def session(self):
        return self.registry.get(self.client.cookies[COOKIE])

    def start(self, field, name, size):
        return self.client.post("/api/uploads", json=dict(field=field, filename=name, size=size))

    def upload(self, field, name, contents=b"uploaded file", chunk_bytes=4, finish=True):
        """The last response: the start's on rejection, else the last chunk's."""
        response = self.start(field, name, len(contents))
        if response.status_code != 200:
            return response
        upload_id = response.json()["upload_id"]
        for index, offset in enumerate(range(0, len(contents), chunk_bytes)):
            if not finish and offset + chunk_bytes >= len(contents):
                break
            response = self.client.put(f"/api/uploads/{field}/{upload_id}/{index}",
                                       content=contents[offset:offset + chunk_bytes])
        return response

    def upload_required(self):
        for field, name in REQUIRED:
            self.assertTrue(self.upload(field, name).json()["complete"])

    def status(self, field):
        return self.client.get("/api/uploads").json()[field]

    def post_load(self, **body):
        body.setdefault("source", "upload")
        if body["source"] == "upload":
            body.setdefault("bins_size", 8)
        return self.client.post("/api/load", json=dict(pseudospots=False, **body))

    def wait_load(self):
        for _ in range(400):
            state = self.client.get("/api/load").json()
            if state.get("state") != "running":
                return state
            time.sleep(.05)
        self.fail("Loading did not settle")

    @contextmanager
    def captured_loads(self):
        """Loads record their LoadConfig and stop (no data is read)."""
        configs = []

        class Recording:
            def __init__(self, config, progress=None, segmentation=None):
                configs.append(config)

            def run(self):
                raise RuntimeError("stopped by the test")

        with patch("hd_viewer.api.routers.load.LoadPipeline", Recording):
            yield configs

    def accepted(self, **body):
        """Whether POST /api/load starts a load (then waits for it)."""
        with self.captured_loads() as configs:
            response = self.post_load(**body)
            if response.status_code == 202:
                self.wait_load()
        self.assertEqual(len(configs), response.status_code == 202, response.text)
        return response.status_code == 202, response.json(), (configs or [None])[0]


class LoadReadinessTests(UploadApiTestCase):
    def test_failed_example_unlocks_settings_and_allows_retry(self):
        with patch("hd_viewer.api.routers.load.LoadPipeline.run", side_effect=ValueError("Example image missing")):
            self.assertEqual(self.post_load(source="example", dataset_id="Tiny").status_code, 202)
            state = self.wait_load()
        self.assertEqual(state["state"], "failed")
        self.assertIn("Example image missing", state["error"])
        self.assertIsNone(self.session.loaded)
        self.assertEqual(self.start("expression", "counts.h5", 4).status_code, 200)   # uploads unlocked
        self.load_example()                                                             # retry works
        self.assertIn("bins", self.session.views)

    def test_viewer_build_failure_also_unlocks_settings(self):
        with patch("hd_viewer.api.routers.layers.open_workspace", side_effect=ValueError("Cannot build viewer")):
            self.assertEqual(self.post_load(source="example", dataset_id="Tiny").status_code, 202)
            state = self.wait_load()
        self.assertEqual(state["state"], "failed")
        self.assertIn("Cannot build viewer", state["error"])
        self.assertEqual(self.start("expression", "counts.h5", 4).status_code, 200)
        self.load_example()

    def test_settings_and_uploads_are_locked_while_loading(self):
        self.upload("expression", "counts.h5", b"counts")
        release, running = threading.Event(), threading.Event()

        def hold():
            running.set()
            release.wait(10)
            raise RuntimeError("released")

        with patch("hd_viewer.api.routers.load.LoadPipeline.run", side_effect=hold):
            self.assertEqual(self.post_load(source="example", dataset_id="Tiny").status_code, 202)
            try:
                self.assertTrue(running.wait(10))
                self.assertEqual(self.post_load(source="example", dataset_id="Tiny").status_code, 409)
                self.assertEqual(self.start("expression", "other.h5", 5).status_code, 409)
                self.assertEqual(self.client.delete("/api/uploads/expression").status_code, 409)
                self.assertEqual(self.session.uploads["expression"].path.read_bytes(), b"counts")
            finally:
                release.set()
            self.assertEqual(self.wait_load()["state"], "failed")
        self.assertEqual(self.client.delete("/api/uploads/expression").status_code, 200)
        self.assertEqual(self.status("expression")["status"], "empty")

    def test_required_files_must_all_arrive(self):
        self.assertFalse(self.accepted()[0])
        self.upload("expression", "counts.h5")
        self.upload("positions", "positions.parquet")
        accepted, reply, _ = self.accepted()
        self.assertFalse(accepted)
        self.assertIn("image", reply["detail"])
        self.upload("image", "tissue.png", b"partially sent", finish=False)
        self.assertFalse(self.accepted()[0])
        self.assertEqual(self.status("image")["status"], "uploading")
        self.assertEqual(self.status("image")["received"], 12)
        self.upload("image", "tissue.png")
        self.assertEqual(self.status("image")["status"], "complete")
        self.assertFalse(self.accepted(bins_size=5)[0])              # the bin size of the files is required
        self.assertTrue(self.accepted()[0])

    def test_failed_load_keeps_readiness_rules(self):
        self.upload_required()
        with patch("hd_viewer.api.routers.load.LoadPipeline.run", side_effect=ValueError("Bad upload")):
            self.assertEqual(self.post_load().status_code, 202)
            self.assertEqual(self.wait_load()["state"], "failed")
        self.assertEqual(self.start("expression", "replacement.h5", 10).status_code, 200)
        accepted, reply, _ = self.accepted()
        self.assertFalse(accepted)
        self.assertIn("expression", reply["detail"])

    def test_optional_uploads_wait_independently(self):
        self.upload_required()
        self.upload("reference", "reference.h5ad", finish=False)
        accepted, reply, _ = self.accepted()
        self.assertFalse(accepted)
        self.assertIn("not finished uploading", reply["detail"])
        self.upload("reference", "reference.h5ad")
        self.assertTrue(self.accepted()[0])

    def test_pending_annotations_block_the_load(self):
        # A complete reference must not load without the annotations still arriving
        self.upload_required()
        self.upload("reference", "reference.h5ad")
        self.upload("reference_annotations", "labels.csv", finish=False)
        self.assertFalse(self.accepted()[0])

    def test_replacement_blocks_load_even_with_previous_contents(self):
        self.upload_required()
        self.upload("reference", "reference.h5ad")
        previous = self.session.uploads["reference"].path
        self.assertTrue(previous.exists())
        self.assertEqual(self.start("reference", "replacement.h5ad", len(b"replacement")).status_code, 200)
        self.assertFalse(previous.exists())
        self.assertFalse(self.accepted()[0])
        self.upload("reference", "replacement.h5ad", b"replacement")
        accepted, _, config = self.accepted()
        self.assertTrue(accepted)
        self.assertEqual(Path(config.reference_path).read_bytes(), b"replacement")

    def test_invalid_extension_is_rejected_and_blocks_load(self):
        self.upload_required()
        response = self.upload("image", "<script>.txt", b"invalid image")
        self.assertEqual(response.status_code, 422)
        self.assertIn("Invalid file", response.json()["detail"])
        self.assertEqual(list(self.session.directory.rglob("*.png")), [])
        self.assertEqual(self.status("image")["status"], "error")
        self.assertIn("Invalid file", self.status("image")["error"])
        self.assertFalse(self.accepted()[0])

    def test_examples_ignore_uploads(self):
        self.upload("expression", "counts.h5", finish=False)
        accepted, _, config = self.accepted(source="example", dataset_id="Tiny")
        self.assertTrue(accepted)
        self.assertNotEqual(config.adata_path, "")
        self.assertFalse(Path(config.adata_path).is_relative_to(settings.UPLOAD_DIR))
        self.assertFalse(self.accepted()[0])                       # your own data still needs every file

    def test_load_reads_the_assembled_files(self):
        contents = {"expression": ("counts.h5", b"expression" * 5),
                    "positions": ("positions.parquet", b"positions" * 3),
                    "image": ("tissue.tif", b"image" * 7)}
        for field, (name, data) in contents.items():
            self.upload(field, name, data, chunk_bytes=6)
        accepted, _, config = self.accepted(bins_size=8)
        self.assertTrue(accepted)
        self.assertEqual(config.bins_size, 8)
        self.assertEqual(config.base_dir, self.session.directory)
        for key, field in (("adata_path", "expression"), ("positions_path", "positions"), ("img_path", "image")):
            name, data = contents[field]
            path = Path(getattr(config, key))
            self.assertTrue(path.is_relative_to(self.session.directory))
            self.assertEqual(path.name, name)
            self.assertEqual(path.read_bytes(), data)
        self.assertIsNone(config.reference_path)


class UploadFieldTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.field = UploadField("expression", (".h5",), lambda: self.folder, max_bytes=100, chunk_bytes=8)

    def send(self, contents, name="counts.h5", upload_id="id"):
        error = self.field.begin(upload_id, name, len(contents))
        for index, offset in enumerate(range(0, len(contents), 8)):
            error = error or self.field.write(index, contents[offset:offset + 8])
        return error

    def test_chunks_are_assembled_in_order_and_renamed_when_complete(self):
        contents = bytes(range(37))
        self.assertIsNone(self.field.begin("id", "counts.h5", len(contents)))
        partial = self.folder / ".expression.part"
        for index, offset in enumerate(range(0, len(contents), 8)):
            self.assertFalse(self.field.complete)
            self.assertTrue(self.field.is_current("id"))
            self.assertIsNone(self.field.write(index, contents[offset:offset + 8]))
            if not self.field.complete:
                self.assertEqual(partial.stat().st_size, offset + 8)
        self.assertTrue(self.field.complete)
        self.assertFalse(partial.exists())
        self.assertEqual(self.field.path, self.folder / "expression" / "counts.h5")
        self.assertEqual(self.field.path.read_bytes(), contents)
        self.assertIn("Uploaded: counts.h5", self.field.describe(True))

    def test_out_of_order_chunk_aborts_the_upload(self):
        self.field.begin("id", "counts.h5", 16)
        self.assertIn("out of order", self.field.write(1, b"12345678"))
        self.assertIn("out of order", self.field.describe(True))
        self.assertEqual(self.field.status, UploadField.ERROR)
        self.assertFalse((self.folder / ".expression.part").exists())
        # Later chunks of the aborted upload are refused.
        self.assertFalse(self.field.is_current("id"))
        self.assertFalse(self.field.complete)

    def test_oversize_files_and_chunks_are_rejected(self):
        self.assertIn("too large", self.field.begin("a", "counts.h5", 101))
        self.assertFalse((self.folder / ".expression.part").exists())

        self.field.begin("b", "counts.h5", 10)
        self.assertIsNone(self.field.write(0, b"12345678"))
        self.assertIn("more data than", self.field.write(1, b"12345678"))

        self.field.begin("c", "counts.h5", 50)
        self.assertIn("invalid chunk size", self.field.write(0, b"x" * 40))
        self.field.begin("d", "counts.h5", 50)
        self.assertIn("invalid chunk size", self.field.write(0, b""))
        self.assertFalse(self.field.complete)

        for size in (0, -1, True, 1.5):
            self.assertIsNotNone(self.field.begin("e", "counts.h5", size))
        with patch("hd_viewer.core.uploads.shutil.disk_usage", return_value=type("U", (), {"free": 5})()):
            self.assertIn("disk space", self.field.begin("f", "counts.h5", 50))
        self.assertIn("Expected .h5", self.field.begin("g", "counts.csv", 5))

    def test_filenames_are_reduced_to_their_basename(self):
        for name in ("../../outside.h5", "/tmp/outside.h5", "..\\..\\outside.h5", "a/b/../outside.h5"):
            with self.subTest(name=name):
                self.assertIsNone(self.send(b"x", name))
                self.assertEqual(self.field.path, self.folder / "expression" / "outside.h5")
                self.assertTrue(self.field.path.resolve().is_relative_to(self.folder.resolve()))
        for name in ("..", "", "/", "bad\x00name.h5"):
            with self.subTest(name=name):
                self.assertIsNotNone(self.send(b"x", name))
                self.assertIsNone(self.field.path)
        self.assertEqual(sorted(p.relative_to(self.folder).as_posix() for p in self.folder.rglob("*")),
                         ["expression"])

    def test_replacing_or_removing_discards_previous_data(self):
        self.send(b"first", "first.h5")
        first = self.field.path
        self.field.begin("second-id", "second.h5", 6)
        self.assertFalse(first.exists())
        self.assertFalse(self.field.complete)
        # A late chunk from an abandoned upload is not accepted for the new one.
        self.assertFalse(self.field.is_current("first-id"))
        self.assertEqual(self.field.status, UploadField.UPLOADING)
        self.assertIsNone(self.field.write(0, b"second"))
        self.assertTrue(self.field.complete)
        self.field.clear()
        self.assertEqual(self.field.status, UploadField.EMPTY)
        self.assertFalse((self.folder / "expression" / "second.h5").exists())

        self.field.begin("third-id", "third.h5", 20)
        self.field.write(0, b"12345678")
        self.field.clear()
        self.assertFalse(self.field.is_current("third-id"))
        self.assertFalse((self.folder / ".expression.part").exists())


class SessionFolderTests(unittest.TestCase):
    def test_client_session_ids_stay_inside_the_upload_folder(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(settings, "UPLOAD_DIR", Path(folder) / "uploads"):
            folders = {session_upload_dir(sid) for sid in ("../../escape", "/absolute/path", "first", "second")}
            self.assertEqual(len(folders), 4)
            for path in folders:
                self.assertEqual(path.parent, Path(folder) / "uploads")


if __name__ == "__main__":
    unittest.main()
