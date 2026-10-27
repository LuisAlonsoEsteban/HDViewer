import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.models import Button, Select

from bokeh_server.panels.startup import StartupController
from src.singlecell import read_reference
from tests.test_deconvolution import ImmediateDoc, tenx_reference_bytes
from tests.test_param_form import make_layer


class ExampleReferenceTests(unittest.TestCase):
    def test_example_reads_local_counts_and_aligns_spatial_metadata_without_uploads(self):
        with tempfile.TemporaryDirectory() as folder:
            reference = Path(folder) / "counts.h5"
            reference.write_bytes(tenx_reference_bytes())
            annotations = Path(folder) / "spatial.parquet"
            pd.DataFrame({"barcode": ["s1", "s0"], "DeconvolutionClass": ["B", "A"]}).to_parquet(annotations)
            cell_annotations = Path(folder) / "SingleCell_MetaData.csv"
            cell_annotations.write_text("Barcode,Patient,Level1,Level2\nc,P1,B,B subtype\nb,P1,A,A subtype\na,P1,A,A subtype\n")
            config = {"reference_path": str(reference), "spatial_annotations_path": str(annotations),
                      "reference_annotations_path": str(cell_annotations), "reference_label_column": "Level1",
                      "image_path": str(Path(folder) / "image.tif"), "expression_path": "unused.h5",
                      "positions_path": "unused.parquet", "bins_size": 8}
            spatial = ad.AnnData(np.ones((2, 2)), obs=pd.DataFrame(index=["s0", "s1"]))
            dataset = SimpleNamespace(adata=spatial, df_tissue_positions=pd.DataFrame(), bins_size=8)
            with patch.dict("bokeh_server.panels.startup.AVAILABLE_DATASETS", {"Human CRC": config}):
                startup = StartupController(ImmediateDoc())
                startup.example_select.value = "Human CRC"
                startup.data_source.active = 1
                self.assertFalse(startup.reference_input_widget.visible)
                self.assertFalse(startup.reference_input.value)
                self.assertFalse(startup.reference_annotations_input.value)
                self.assertIsNone(startup._reference_upload_error())
                with patch("bokeh_server.panels.startup.DataSet", return_value=dataset), \
                     patch.object(startup, "generate_tiles"), \
                     patch.object(startup, "generate_pseudospots"), \
                     patch.object(startup, "generate_segmentations"), \
                     patch.object(startup, "build_b2s") as build, \
                     patch("src.singlecell.NamedTemporaryFile", side_effect=AssertionError("Local files must not be copied")):
                    startup._prepare_data()
                self.assertEqual(startup.app.deconvolution_reference.n_obs, 3)
                self.assertEqual(startup.app.deconvolution_reference_name, "counts.h5")
                self.assertEqual(startup.app.deconvolution_reference.obs.Level1.tolist(), ["A", "A", "B"])
                layer, _ = make_layer()
                startup.app.add_geometry(layer)
                widgets = layer.deconvolution_manager.get_widgets(ImmediateDoc())
                labels = next(w for w in widgets if isinstance(w, Select))
                self.assertEqual(labels.options, ["Patient", "Level1", "Level2"])
                self.assertEqual(labels.value, "Level1")
                self.assertFalse(labels.disabled)
                self.assertFalse(next(w for w in widgets if isinstance(w, Button)).disabled)
                self.assertEqual(spatial.obs.DeconvolutionClass.tolist(), ["A", "B"])
                build.assert_called_once_with(spatial, dataset.df_tissue_positions, 8)
                # Every example hides uploads; your own data can still supply a reference.
                startup.example_select.value = "Mouse brain (8um)"
                self.assertFalse(startup.reference_input_widget.visible)
                startup.data_source.active = 0
                self.assertTrue(startup.reference_input_widget.visible)

    def test_example_without_reference_ignores_stale_uploads_and_loads_normally(self):
        startup = StartupController(ImmediateDoc())
        startup.reference_input.set_from_json("filename", "pending.h5")
        startup.reference_annotations_input.set_from_json("value", "stale annotation upload")
        startup.app.deconvolution_reference = object()
        startup.app.deconvolution_reference_name = "previous.h5"
        startup.app.deconvolution_label_column = "old_label"
        startup.data_source.active = 1
        startup.example_select.value = "Mouse brain (8um)"
        self.assertFalse(startup.reference_input_widget.visible)
        self.assertIsNone(startup._reference_upload_error())
        with patch("src.singlecell.read_reference") as read, \
             patch("bokeh_server.panels.startup.DataSet"), \
             patch.object(startup, "generate_tiles") as tiles, \
             patch.object(startup, "generate_pseudospots"), \
             patch.object(startup, "generate_segmentations"), \
             patch.object(startup, "ensure_b2s") as bins:
            startup._prepare_data()
        read.assert_not_called()
        tiles.assert_called_once()
        bins.assert_called_once()
        self.assertIsNone(startup.app.deconvolution_reference)
        self.assertIsNone(startup.app.deconvolution_reference_name)
        self.assertEqual(startup.app.deconvolution_label_column, "")
        layer, _ = make_layer()
        startup.app.add_geometry(layer)
        widgets = layer.deconvolution_manager.get_widgets(ImmediateDoc())
        self.assertTrue(next(w for w in widgets if isinstance(w, Button)).disabled)
        self.assertTrue(any("No single-cell reference" in getattr(w, "text", "") for w in widgets))

    def test_local_reference_annotation_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            reference = Path(folder) / "counts.h5"
            reference.write_bytes(tenx_reference_bytes())
            annotations = Path(folder) / "cells.parquet"
            pd.DataFrame({"barcode": ["c", "b", "a"], "cell_type": ["B", "A", "A"]}).to_parquet(annotations)
            with patch("src.singlecell.NamedTemporaryFile", side_effect=AssertionError("Unexpected file copy")):
                result = read_reference(reference, annotations)
            self.assertEqual(result.obs.cell_type.tolist(), ["A", "A", "B"])


if __name__ == "__main__":
    unittest.main()
