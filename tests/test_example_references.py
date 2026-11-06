"""Example datasets read their single-cell reference and annotations from local paths."""
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from hd_viewer.startup import datasets
from hd_viewer.startup.data.singlecell import read_reference
from tests.test_deconvolution import DeconvolutionApiTestCase, tenx_reference_bytes


class ExampleReferenceTests(DeconvolutionApiTestCase):
    def test_example_reads_local_counts_and_aligns_spatial_metadata_without_uploads(self):
        folder = Path(self.folder.name)
        reference = folder / "counts.h5"
        reference.write_bytes(tenx_reference_bytes())
        annotations = folder / "spatial.parquet"
        # Spatial barcodes of the Tiny dataset are a, b, c: aligned by barcode, not row order
        pd.DataFrame({"barcode": ["c", "a", "b"], "DeconvolutionClass": ["C", "A", "B"]}).to_parquet(annotations)
        cell_annotations = folder / "SingleCell_MetaData.csv"
        cell_annotations.write_text("Barcode,Patient,Level1,Level2\nc,P1,B,B subtype\nb,P1,A,A subtype\n"
                                    "a,P1,A,A subtype\n")
        datasets.AVAILABLE_DATASETS["Tiny"] = dict(
            datasets.AVAILABLE_DATASETS["Tiny"], reference_path=str(reference),
            spatial_annotations_path=str(annotations), reference_annotations_path=str(cell_annotations),
            reference_label_column="Level1")
        with patch("hd_viewer.startup.data.singlecell.NamedTemporaryFile",
                   side_effect=AssertionError("Local files must not be copied")):
            self.load_example()
        self.assertEqual(self.client.get("/api/uploads").json()["reference"]["status"], "empty")
        loaded = self.session.loaded
        self.assertEqual(loaded.reference.n_obs, 3)
        self.assertEqual(loaded.reference_name, "counts.h5")
        self.assertEqual(loaded.reference.obs.Level1.tolist(), ["A", "A", "B"])
        obs = self.session.views["bins"].sample.adata.obs
        self.assertEqual(obs.DeconvolutionClass.astype(str).to_dict(), {"a": "A", "b": "B", "c": "C"})

        state = self.state()
        self.assertEqual(state["unavailable"], "")
        self.assertEqual(state["columns"], ["Patient", "Level1", "Level2"])
        self.assertEqual(state["label_column"], "Level1")
        self.assertEqual(state["reference"]["name"], "counts.h5")
        self.assertEqual(self.client.post("/api/layers/bins/deconvolution/check",
                                          json=dict(method="nnls", label_column="Level1")).status_code, 200)

    def test_example_without_reference_ignores_stale_uploads_and_loads_normally(self):
        self.upload("reference", "pending.h5", finish=False)
        self.upload("reference_annotations", "labels.csv", b"stale annotation upload")
        with patch("hd_viewer.startup.data.singlecell.read_reference") as read:
            self.load_example()
        read.assert_not_called()
        loaded = self.session.loaded
        self.assertIsNone(loaded.reference)
        self.assertIsNone(loaded.reference_name)
        self.assertEqual(loaded.reference_label_column, "")
        self.assertIn("bins", self.session.views)
        state = self.state()
        self.assertIn("No single-cell reference", state["unavailable"])
        self.assertIsNone(state["reference"])
        self.assertEqual(self.client.post("/api/layers/bins/deconvolution",
                                          json=dict(method="nnls", label_column="cell_type")).status_code, 422)

    def test_local_reference_annotation_paths(self):
        folder = Path(self.folder.name)
        reference = folder / "counts.h5"
        reference.write_bytes(tenx_reference_bytes())
        annotations = folder / "cells.parquet"
        pd.DataFrame({"barcode": ["c", "b", "a"], "cell_type": ["B", "A", "A"]}).to_parquet(annotations)
        with patch("hd_viewer.startup.data.singlecell.NamedTemporaryFile",
                   side_effect=AssertionError("Unexpected file copy")):
            result = read_reference(reference, annotations)
        self.assertEqual(result.obs.cell_type.tolist(), ["A", "A", "B"])


if __name__ == "__main__":
    unittest.main()
