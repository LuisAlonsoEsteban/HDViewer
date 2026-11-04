import base64
import json
import tempfile
import unittest
import threading
import time
import uuid
from pathlib import Path
from unittest.mock import patch

from bokeh.document import Document
from bokeh.models import Div
from bokeh.server.callbacks import NextTickCallback

from hd_viewer.startup import StartupController
from hd_viewer.startup.uploads import ChunkedUpload


def send(field, **message):
    """Apply a browser message to the field's channel, as a websocket patch would."""
    field.channel.set_from_json("request", json.dumps(message))
    return json.loads(field.channel.response)


def start_upload(field, name, size):
    upload_id = uuid.uuid4().hex
    return upload_id, send(field, action="start", id=upload_id, name=name, size=size)


def send_chunk(field, upload_id, index, data):
    return send(field, action="chunk", id=upload_id, index=index,
                data=base64.b64encode(data).decode())


def upload_file(field, name, contents=b"uploaded file", chunk_bytes=None, finish=True):
    """Drive the chunked-upload protocol exactly like static/js/chunked_upload.js."""
    upload_id, reply = start_upload(field, name, len(contents))
    if not reply["ok"]:
        return reply
    chunk_bytes = chunk_bytes or field.chunk_bytes
    for index, offset in enumerate(range(0, len(contents), chunk_bytes)):
        if not finish and offset + chunk_bytes >= len(contents):
            break
        reply = send_chunk(field, upload_id, index, contents[offset:offset + chunk_bytes])
    return reply


def use_temporary_upload_dir(test, startup):
    folder = tempfile.TemporaryDirectory()
    test.addCleanup(folder.cleanup)
    startup.SESSION_UPLOAD_DIR = Path(folder.name) / "session"
    return startup.SESSION_UPLOAD_DIR


class StartupUploadTests(unittest.TestCase):
    def setUp(self):
        def skip_segmentation(controller):
            controller.segmentation_controls = {}
            controller.segmentation_info = Div()

        with patch.object(StartupController, "_build_segmentation_widgets", skip_segmentation), \
                patch("hd_viewer.startup.App"):
            self.startup = StartupController(Document())
        self.upload_dir = use_temporary_upload_dir(self, self.startup)

    def receive(self, field, filename, contents=b"uploaded file"):
        reply = upload_file(field, filename, contents, chunk_bytes=4)
        self.assertEqual(reply["done"], reply["ok"])
        return reply

    def receive_required(self):
        self.receive(self.startup.h5_input, "counts.h5")
        self.receive(self.startup.positions_input, "positions.parquet")
        self.receive(self.startup.image_input, "tissue.png")

    def settle_callbacks(self, predicate):
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            for callback in list(self.startup.doc.session_callbacks):
                if isinstance(callback, NextTickCallback):
                    self.startup.doc.remove_next_tick_callback(callback)
                    callback.callback()
            if predicate():
                return
            time.sleep(.005)
        self.fail("Loading did not settle")

    def test_failed_example_unlocks_settings_and_allows_retry(self):
        self.startup.data_source.active = 1
        old_data, old_bins, old_image = object(), object(), object()

        def fail():
            self.assertTrue(all(widget.disabled for widget in self.startup._loading_controls))
            self.startup.dataset = old_data
            self.startup.b2s = old_bins
            self.startup.image_processor = old_image
            raise ValueError("Example image missing")

        with patch.object(self.startup, '_prepare_data', side_effect=fail):
            self.startup.load_data_callback()
            self.settle_callbacks(lambda: not self.startup.loading)
        self.assertTrue(self.startup.failed)
        self.assertTrue(all(not widget.disabled for widget in self.startup._loading_controls))
        self.assertIsNone(self.startup._progress_timer)
        self.assertIn('Example image missing', self.startup.status_text.text)
        self.assertIn('try again', self.startup.status_text.text)
        self.startup.pseudospots_checkbox.active = False

        def succeed():
            self.assertFalse(self.startup.failed)
            self.assertFalse(self.startup.finished)
            self.assertIsNone(self.startup.dataset)
            self.assertIsNone(self.startup.b2s)
            self.assertIsNone(self.startup.image_processor)
            for task in self.startup.tasks:
                self.startup.tasks[task] = True

        with patch.object(self.startup, '_prepare_data', side_effect=succeed):
            self.startup.load_data_callback()
            self.settle_callbacks(lambda: self.startup._ready_to_finish)
        self.startup._finish_loading()
        self.assertTrue(self.startup.finished)
        self.assertEqual(self.startup.status_text.text, 'Done')
        self.startup.app.build.assert_called_once()

    def test_viewer_build_failure_also_unlocks_settings(self):
        self.startup.data_source.active = 1
        self.startup.app.build.side_effect = ValueError('Cannot build viewer')

        def prepare():
            for task in self.startup.tasks:
                self.startup.tasks[task] = True

        with patch.object(self.startup, '_prepare_data', side_effect=prepare):
            self.startup.load_data_callback()
            self.settle_callbacks(lambda: not self.startup.loading)
        self.assertFalse(self.startup.finished)
        self.assertFalse(self.startup.load_button.disabled)
        self.assertFalse(self.startup.example_select.disabled)
        self.assertIn('Cannot build viewer', self.startup.status_text.text)

    def test_failed_upload_keeps_load_readiness_rules(self):
        self.receive_required()
        with patch.object(self.startup, '_prepare_data', side_effect=ValueError('Bad upload')):
            self.startup.load_data_callback()
            self.settle_callbacks(lambda: not self.startup.loading)
        self.assertFalse(self.startup.h5_input.disabled)
        self.assertFalse(self.startup.load_button.disabled)
        self.startup.data_source.active = 0
        start_upload(self.startup.h5_input, 'replacement.h5', 10)
        self.assertTrue(self.startup.load_button.disabled)

    def test_failure_waits_for_other_workers_before_enabling_retry(self):
        self.startup.data_source.active = 1
        release = threading.Event()

        def prepare():
            self.startup.run_in_thread(lambda: release.wait(5), task='spots')
            raise ValueError('Preparation failed')

        try:
            with patch.object(self.startup, '_prepare_data', side_effect=prepare):
                self.startup.load_data_callback()
                self.settle_callbacks(lambda: self.startup.failed and not self.startup.preparing)
            self.assertTrue(self.startup.load_button.disabled)
            self.assertTrue(self.startup.loading)
        finally:
            release.set()
        self.settle_callbacks(lambda: not self.startup.loading)
        self.assertFalse(self.startup.load_button.disabled)
        self.assertEqual(self.startup._running_workers, 0)

    def test_required_files_must_all_arrive(self):
        self.assertTrue(self.startup.load_button.disabled)
        self.receive(self.startup.h5_input, "counts.h5")
        self.receive(self.startup.positions_input, "positions.parquet")
        self.assertTrue(self.startup.load_button.disabled)
        upload_file(self.startup.image_input, "tissue.png", b"partially sent", chunk_bytes=4, finish=False)
        self.assertTrue(self.startup.load_button.disabled)
        self.assertIn("Uploading tissue.png", self.startup.upload_status.text)
        self.receive(self.startup.image_input, "tissue.png")
        self.assertFalse(self.startup.load_button.disabled)
        self.assertIn("Uploaded: tissue.png", self.startup.upload_status.text)

    def test_optional_uploads_wait_independently(self):
        self.receive_required()
        for field, name in [(self.startup.reference_input, "reference.h5ad"),
                            (self.startup.reference_annotations_input, "labels.csv")]:
            start_upload(field, name, 10)
        self.assertTrue(self.startup.load_button.disabled)
        self.receive(self.startup.reference_input, "reference.h5ad")
        self.assertTrue(self.startup.load_button.disabled)
        self.assertIn("Uploading labels.csv", self.startup.upload_status.text)
        self.receive(self.startup.reference_annotations_input, "labels.csv")
        self.assertFalse(self.startup.load_button.disabled)

    def test_replacement_blocks_load_even_with_previous_contents(self):
        self.receive_required()
        self.receive(self.startup.reference_input, "reference.h5ad")
        start_upload(self.startup.reference_input, "replacement.h5ad", len(b"replacement"))
        self.assertFalse((self.upload_dir / "reference" / "reference.h5ad").exists())
        self.assertTrue(self.startup.load_button.disabled)
        with patch.object(self.startup, "restart") as restart:
            self.startup.load_data_callback()
        restart.assert_not_called()
        self.receive(self.startup.reference_input, "replacement.h5ad", b"replacement")
        self.assertFalse(self.startup.load_button.disabled)

    def test_invalid_extension_blocks_load_and_escapes_filename(self):
        self.receive_required()
        reply = upload_file(self.startup.image_input, "<script>.txt", b"invalid image")
        self.assertFalse(reply["ok"])
        self.assertEqual(list(self.upload_dir.rglob("*.png")), [])
        self.assertTrue(self.startup.load_button.disabled)
        self.assertIn("Invalid file", self.startup.upload_status.text)
        self.assertIn("&lt;script&gt;", self.startup.upload_status.text)

    def test_examples_ignore_uploads(self):
        self.startup.data_source.active = 1
        self.assertFalse(self.startup.load_button.disabled)
        start_upload(self.startup.h5_input, "counts.h5", 10)
        self.assertFalse(self.startup.load_button.disabled)
        self.startup.data_source.active = 0
        self.assertTrue(self.startup.load_button.disabled)

class ChunkedUploadTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.changes = 0
        self.locked = False
        self.field = ChunkedUpload("expression", (".h5",), lambda: self.folder,
                                   on_change=self.count_change, can_upload=lambda: not self.locked,
                                   max_bytes=100, chunk_bytes=8)

    def count_change(self):
        self.changes += 1

    def test_chunks_are_assembled_in_order_and_renamed_when_complete(self):
        contents = bytes(range(37))
        upload_id, reply = start_upload(self.field, "counts.h5", len(contents))
        self.assertTrue(reply["ok"])
        partial = self.folder / ".expression.part"
        for index, offset in enumerate(range(0, len(contents), 8)):
            self.assertFalse(self.field.complete)
            reply = send_chunk(self.field, upload_id, index, contents[offset:offset + 8])
            self.assertTrue(reply["ok"])
            self.assertEqual(reply["index"], index)
            # The chunk is dropped from the model once written.
            self.assertEqual(self.field.channel.request, "")
            self.assertEqual(partial.stat().st_size if not reply["done"] else len(contents),
                             min(offset + 8, len(contents)))
        self.assertTrue(reply["done"])
        self.assertTrue(self.field.complete)
        self.assertFalse(partial.exists())
        self.assertEqual(self.field.path, self.folder / "expression" / "counts.h5")
        self.assertEqual(self.field.path.read_bytes(), contents)
        self.assertIn("Uploaded: counts.h5", self.field.describe(True))
        self.assertTrue(self.field.remove_button.visible)

    def test_out_of_order_chunk_aborts_the_upload(self):
        upload_id, _ = start_upload(self.field, "counts.h5", 16)
        reply = send_chunk(self.field, upload_id, 1, b"12345678")
        self.assertFalse(reply["ok"])
        self.assertIn("out of order", self.field.describe(True))
        self.assertEqual(self.field.status, ChunkedUpload.ERROR)
        self.assertFalse((self.folder / ".expression.part").exists())
        # Later chunks of the aborted upload are ignored.
        self.assertFalse(send_chunk(self.field, upload_id, 0, b"12345678")["ok"])
        self.assertFalse(self.field.complete)

    def test_oversize_files_and_chunks_are_rejected(self):
        upload_id, reply = start_upload(self.field, "counts.h5", 101)
        self.assertFalse(reply["ok"])
        self.assertIn("too large", reply["error"])
        self.assertFalse((self.folder / ".expression.part").exists())

        upload_id, _ = start_upload(self.field, "counts.h5", 10)
        reply = send_chunk(self.field, upload_id, 0, b"12345678")
        reply = send_chunk(self.field, upload_id, 1, b"12345678")
        self.assertFalse(reply["ok"])
        self.assertIn("more data than", reply["error"])

        upload_id, _ = start_upload(self.field, "counts.h5", 50)
        reply = send_chunk(self.field, upload_id, 0, b"x" * 40)
        self.assertFalse(reply["ok"])
        self.assertIn("invalid chunk size", reply["error"])

        upload_id, _ = start_upload(self.field, "counts.h5", 50)
        reply = send(self.field, action="chunk", id=upload_id, index=0, data="not base64!")
        self.assertFalse(reply["ok"])
        self.assertFalse(self.field.complete)

        upload_id, reply = start_upload(self.field, "counts.h5", 0)
        self.assertFalse(reply["ok"])

    def test_filenames_are_reduced_to_their_basename(self):
        for name in ("../../outside.h5", "/tmp/outside.h5", "..\\..\\outside.h5", "a/b/../outside.h5"):
            with self.subTest(name=name):
                reply = upload_file(self.field, name, b"x")
                self.assertTrue(reply["done"])
                self.assertEqual(self.field.path, self.folder / "expression" / "outside.h5")
                self.assertTrue(self.field.path.resolve().is_relative_to(self.folder.resolve()))
        for name in ("..", "", "/", "bad\x00name.h5"):
            with self.subTest(name=name):
                self.assertFalse(upload_file(self.field, name, b"x")["ok"])
                self.assertIsNone(self.field.path)
        self.assertEqual(sorted(p.relative_to(self.folder).as_posix() for p in self.folder.rglob("*")),
                         ["expression"])

    def test_replacing_or_removing_discards_previous_data(self):
        upload_file(self.field, "first.h5", b"first")
        first = self.field.path
        upload_id, _ = start_upload(self.field, "second.h5", 6)
        self.assertFalse(first.exists())
        self.assertFalse(self.field.complete)
        # A late chunk from an abandoned upload does not touch the new one.
        self.assertFalse(send_chunk(self.field, "old-id", 0, b"zz")["ok"])
        self.assertEqual(self.field.status, ChunkedUpload.UPLOADING)
        self.assertTrue(send_chunk(self.field, upload_id, 0, b"second")["done"])
        self.field.remove()
        self.assertEqual(self.field.status, ChunkedUpload.EMPTY)
        self.assertFalse((self.folder / "expression" / "second.h5").exists())
        self.assertFalse(self.field.remove_button.visible)

        upload_id, _ = start_upload(self.field, "third.h5", 20)
        send_chunk(self.field, upload_id, 0, b"12345678")
        self.field.remove()
        response = json.loads(self.field.channel.response)
        self.assertEqual((response["id"], response["ok"]), (upload_id, False))
        self.assertFalse((self.folder / ".expression.part").exists())

    def test_uploads_are_locked_while_loading(self):
        upload_file(self.field, "counts.h5", b"counts")
        self.locked = True
        reply = upload_file(self.field, "other.h5", b"other")
        self.assertFalse(reply["ok"])
        self.field.remove()
        self.assertTrue(self.field.complete)
        self.assertEqual(self.field.path.read_bytes(), b"counts")


class StartupLoadsAssembledFilesTests(unittest.TestCase):
    def test_load_reads_the_assembled_files(self):
        def skip_segmentation(controller):
            controller.segmentation_controls = {}
            controller.segmentation_info = Div()

        with patch.object(StartupController, "_build_segmentation_widgets", skip_segmentation), \
                patch("hd_viewer.startup.App"):
            startup = StartupController(Document())
        upload_dir = use_temporary_upload_dir(self, startup)
        contents = {startup.h5_input: ("counts.h5", b"expression" * 5),
                    startup.positions_input: ("positions.parquet", b"positions" * 3),
                    startup.image_input: ("tissue.tif", b"image" * 7)}
        for field, (name, data) in contents.items():
            upload_file(field, name, data, chunk_bytes=6)
        startup.resolution_select.value = "8"
        self.assertFalse(startup.load_button.disabled)
        with patch("hd_viewer.startup.VisiumHD") as dataset, \
                patch.object(startup, "generate_tiles"), \
                patch.object(startup, "generate_pseudospots"), \
                patch.object(startup, "generate_segmentations"), \
                patch.object(startup, "ensure_b2s"):
            startup._prepare_data()
        kwargs = dataset.call_args.kwargs
        self.assertEqual(kwargs["bins_size"], 8)
        for key, field in (("adata_path", startup.h5_input), ("positions_path", startup.positions_input),
                           ("img_path", startup.image_input)):
            name, data = contents[field]
            path = Path(kwargs[key])
            self.assertTrue(path.is_relative_to(upload_dir))
            self.assertEqual(path.name, name)
            self.assertEqual(path.read_bytes(), data)


if __name__ == "__main__":
    unittest.main()
