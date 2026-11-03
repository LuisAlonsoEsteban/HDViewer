import base64
import unittest
import threading
import time
from unittest.mock import patch

from bokeh.document import Document
from bokeh.models import Div
from bokeh.server.callbacks import NextTickCallback

from hd_viewer.startup import StartupController


class StartupUploadTests(unittest.TestCase):
    def setUp(self):
        def skip_segmentation(controller):
            controller.segmentation_controls = {}
            controller.segmentation_info = Div()

        with patch.object(StartupController, "_build_segmentation_widgets", skip_segmentation), \
                patch("hd_viewer.startup.App"):
            self.startup = StartupController(Document())

    def receive(self, widget, filename, contents=b"uploaded file"):
        widget.set_from_json("filename", filename)
        widget.set_from_json("value", base64.b64encode(contents).decode())

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
        self.startup.data_source.active = 2
        self.assertTrue(self.startup.load_button.disabled)
        self.startup.data_source.active = 0
        self.startup._upload_pending[self.startup.h5_input.id].value = 'replacement.h5'
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
        self.startup.image_input.set_from_json("filename", "tissue.png")
        self.assertTrue(self.startup.load_button.disabled)
        self.receive(self.startup.image_input, "tissue.png")
        self.assertFalse(self.startup.load_button.disabled)
        self.assertIn("Uploaded: tissue.png", self.startup.upload_status.text)

    def test_optional_uploads_wait_independently(self):
        self.receive_required()
        for widget, name in [(self.startup.reference_input, "reference.h5ad"),
                             (self.startup.reference_annotations_input, "labels.csv")]:
            self.startup._upload_pending[widget.id].value = name
        self.assertTrue(self.startup.load_button.disabled)
        self.receive(self.startup.reference_input, "reference.h5ad")
        self.assertTrue(self.startup.load_button.disabled)
        self.assertIn("Uploading labels.csv", self.startup.upload_status.text)
        self.receive(self.startup.reference_annotations_input, "labels.csv")
        self.assertFalse(self.startup.load_button.disabled)

    def test_replacement_blocks_load_even_with_previous_contents(self):
        self.receive_required()
        self.receive(self.startup.reference_input, "reference.h5ad")
        self.startup._upload_pending[self.startup.reference_input.id].value = "replacement.h5ad"
        self.assertTrue(self.startup.load_button.disabled)
        with patch.object(self.startup, "restart") as restart:
            self.startup.load_data_callback()
        restart.assert_not_called()
        self.receive(self.startup.reference_input, "replacement.h5ad", b"replacement")
        self.assertFalse(self.startup.load_button.disabled)

    def test_invalid_extension_blocks_load_and_escapes_filename(self):
        self.receive_required()
        self.receive(self.startup.image_input, "<script>.txt", b"invalid image")
        self.assertTrue(self.startup.load_button.disabled)
        self.assertIn("Invalid file", self.startup.upload_status.text)
        self.assertIn("&lt;script&gt;", self.startup.upload_status.text)

    def test_examples_ignore_uploads_and_url_stays_disabled(self):
        self.startup.data_source.active = 1
        self.assertFalse(self.startup.load_button.disabled)
        self.startup._upload_pending[self.startup.h5_input.id].value = "counts.h5"
        self.assertFalse(self.startup.load_button.disabled)
        self.startup.data_source.active = 2
        self.assertTrue(self.startup.load_button.disabled)
        self.startup.data_source.active = 0
        self.assertTrue(self.startup.load_button.disabled)


if __name__ == "__main__":
    unittest.main()
