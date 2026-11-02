import base64
import unittest
from unittest.mock import patch

from bokeh.document import Document
from bokeh.models import Div

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
