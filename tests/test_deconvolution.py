import tempfile
import base64
import time
import unittest
from io import BytesIO
from unittest.mock import patch
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import Button, FileInput, Select, TextInput
from bokeh.plotting import figure
from scipy import sparse

from hd_viewer.board.operations.layer.deconvolution import deconvolve
from hd_viewer.startup.data.singlecell import read_reference
from hd_viewer.startup.sample import Cells, Spots
from hd_viewer.board.components.layer import Circles, Patches
from hd_viewer.startup import StartupController
from tests.test_param_form import make_layer


def reference():
    # Deliberately reversed genes and unequal library sizes.
    return ad.AnnData(sparse.csr_matrix([[0., 100.], [0., 10.], [40., 0.]]),
                      obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["a", "b", "c"]),
                      var=pd.DataFrame(index=["g1", "g0"]))


def tenx_reference_bytes():
    ref = reference()
    matrix = ref.X.T.tocsc()
    contents = BytesIO()
    with h5py.File(contents, "w") as file:
        group = file.create_group("matrix")
        for name, data in (("data", matrix.data), ("indices", matrix.indices),
                           ("indptr", matrix.indptr), ("shape", matrix.shape),
                           ("barcodes", np.asarray(ref.obs_names, dtype="S"))):
            group.create_dataset(name, data=data)
        features = group.create_group("features")
        for name, data in (("name", np.asarray(ref.var_names, dtype="S")),
                           ("id", np.asarray(ref.var_names, dtype="S")),
                           ("feature_type", [b"Gene Expression"] * ref.n_vars)):
            features.create_dataset(name, data=data)
    return contents.getvalue()


class ImmediateDoc(Document):
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


class DeconvolutionTests(unittest.TestCase):
    def test_pending_reference_upload_keeps_startup_form_usable(self):
        startup = StartupController(ImmediateDoc())
        startup.reference_annotations_input = FileInput(
            value=base64.b64encode(b"barcode,cell_type\na,A\n").decode(), filename="labels.csv")
        with patch.object(startup, "restart") as restart, patch.object(startup, "_prepare_data") as prepare:
            startup.load_data_callback()
            restart.assert_not_called()
            prepare.assert_not_called()
        self.assertFalse(startup.loading)
        self.assertFalse(startup.load_button.disabled)
        self.assertFalse(startup.reference_input.disabled)
        self.assertIn("has not reached the server", startup.status_text.text)

        # Simulate the browser finishing and synchronising the delayed reference upload.
        startup.reference_input.set_from_json("filename", "reference.h5")
        startup.reference_input.set_from_json("value", base64.b64encode(tenx_reference_bytes()).decode())
        self.assertIn("Reference uploaded: reference.h5", startup.reference_upload_status.text)
        self.assertIsNone(startup._reference_upload_error())

    def test_selected_reference_without_contents_is_not_silently_skipped(self):
        startup = StartupController(ImmediateDoc())
        startup.reference_input.set_from_json("filename", "reference.h5")
        with self.assertRaisesRegex(ValueError, "has not reached the server"):
            startup._prepare_deconvolution_reference()

    def test_startup_reference_is_read_once_and_reused_by_later_layers(self):
        startup = StartupController(ImmediateDoc())
        self.assertTrue(startup.reference_input_widget.visible)
        self.assertIn(startup.reference_input, startup.layout.references())
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "reference.h5"
            uploaded_reference = reference()
            uploaded_reference.obs["donor"] = ["d1", "d1", "d2"]
            uploaded_reference.write_h5ad(path)
            startup.reference_input = FileInput(
                value=base64.b64encode(path.read_bytes()).decode(), filename=path.name)
        # Exercise the actual preparation entry point without loading spatial data.
        with patch("anndata.read_h5ad", wraps=ad.read_h5ad) as read, \
             patch.object(startup, "_save_uploaded_file", return_value="unused"), \
             patch("hd_viewer.startup.VisiumHD"), \
             patch.object(startup, "generate_tiles"), \
             patch.object(startup, "generate_pseudospots"), \
             patch.object(startup, "generate_segmentations"), \
             patch.object(startup, "ensure_b2s"):
            startup._prepare_data()
            shared = startup.app.deconvolution_reference
            self.assertIsNotNone(shared)
            for _ in range(2):
                layer, _ = make_layer()
                startup.app.add_geometry(layer)
                self.assertIs(layer.deconvolution_manager.reference, shared)
                tabs = layer.get_widgets(ImmediateDoc())
                widgets = next(t.child.children for t in tabs.tabs if t.title == "Deconvolution")
                run = next(w for w in widgets if isinstance(w, Button) and w.label == "Run deconvolution")
                labels = next(w for w in widgets if isinstance(w, Select))
                self.assertFalse(run.disabled)
                self.assertEqual(labels.value, "cell_type")
                self.assertEqual(labels.options, ["cell_type", "donor"])
                self.assertFalse(labels.disabled)
                self.assertFalse(any(isinstance(w, TextInput) for w in widgets))
                self.assertEqual([w.label for w in widgets if isinstance(w, Button)], ["Run deconvolution"])
                self.assertTrue(any("Reference ready: reference.h5" in getattr(w, "text", "")
                                    for w in widgets))
                run._trigger_event(ButtonClick(run))
                deadline = time.monotonic() + 10
                while run.disabled and time.monotonic() < deadline:
                    time.sleep(.01)
                self.assertFalse(run.disabled)
                self.assertIn("deconvolution", layer.sample.adata.obsm)
                self.assertEqual(layer.sample.adata.uns["deconvolution"]["reference_name"], "reference.h5")
            self.assertEqual(read.call_count, 1)

    def test_tenx_h5_annotations_populate_columns_at_startup(self):
        # Labels must align by barcode, not by annotation row order.
        metadata = pd.DataFrame({"barcode": ["c", "a", "b"], "cell_type": ["B", "A", "A"],
                                 "donor": ["d2", "d1", "d1"]})
        formats = [("csv", metadata.to_csv(index=False).encode())]
        for layout, frame in (("column", metadata), ("index", metadata.set_index("barcode")),
                              ("first_column", metadata.rename(columns={"barcode": "cell_id"}))):
            buffer = BytesIO()
            frame.to_parquet(buffer)
            formats.append((f"{layout}.parquet", buffer.getvalue()))
        for format_name, contents in formats:
            with self.subTest(format=format_name):
                startup = StartupController(ImmediateDoc())
                self.assertIn(".parquet", startup.reference_annotations_input.accept)
                startup.reference_input = FileInput(
                    value=base64.b64encode(tenx_reference_bytes()).decode(), filename="counts.h5")
                startup.reference_annotations_input = FileInput(
                    value=base64.b64encode(contents).decode(), filename=f"annotations.{format_name}")
                startup._prepare_deconvolution_reference()
                layer, _ = make_layer()
                startup.app.add_geometry(layer)
                widgets = layer.deconvolution_manager.get_widgets(ImmediateDoc())
                labels = next(w for w in widgets if isinstance(w, Select))
                run = next(w for w in widgets if isinstance(w, Button))
                self.assertEqual(labels.options, ["cell_type", "donor"])
                self.assertFalse(labels.disabled)
                self.assertFalse(run.disabled)
                self.assertEqual(startup.app.deconvolution_reference.obs.cell_type.tolist(), ["A", "A", "B"])
                np.testing.assert_allclose(startup.app.deconvolution_reference.X.toarray(), reference().X.toarray())

    def test_parquet_annotations_validate_barcodes(self):
        for barcodes, message in ((["a"], "missing 2 reference cell barcodes"),
                                   (["a", "a"], "unique cell barcodes")):
            with self.subTest(barcodes=barcodes):
                buffer = BytesIO()
                pd.DataFrame({"barcode": barcodes, "cell_type": ["A"] * len(barcodes)}).to_parquet(buffer)
                with self.assertRaisesRegex(ValueError, message):
                    read_reference(tenx_reference_bytes(), buffer.getvalue())

    def test_counts_only_h5_explains_missing_annotations(self):
        layer, _ = make_layer()
        layer.deconvolution_manager.reference = read_reference(tenx_reference_bytes())
        widgets = layer.deconvolution_manager.get_widgets(ImmediateDoc())
        self.assertFalse(any(isinstance(w, TextInput) for w in widgets))
        self.assertTrue(next(w for w in widgets if isinstance(w, Button)).disabled)
        self.assertTrue(any("no .obs annotation columns" in getattr(w, "text", "") for w in widgets))
        with self.assertRaisesRegex(ValueError, "missing 2 reference cell barcodes"):
            read_reference(tenx_reference_bytes(), b"barcode,cell_type\na,A\n")
        with self.assertRaisesRegex(ValueError, "unique cell barcodes"):
            read_reference(tenx_reference_bytes(), b"barcode,cell_type\na,A\na,B\n")

    def test_optional_or_invalid_startup_reference(self):
        startup = StartupController(ImmediateDoc())
        startup.app.deconvolution_reference = reference()
        startup._prepare_deconvolution_reference()
        self.assertIsNone(startup.app.deconvolution_reference)
        startup.reference_input = FileInput(
            value=base64.b64encode(b"not an h5ad file").decode(), filename="bad.h5ad")
        with patch.object(startup, "generate_tiles") as tiles:
            with self.assertRaisesRegex(ValueError, "Cannot load single-cell reference bad.h5ad"):
                startup._prepare_data()
            tiles.assert_not_called()

    def test_tab_is_available_on_spots_bins_and_cell_methods(self):
        bins, _ = make_layer()
        for component_class, sample_class, alias in ((Circles, Spots, "Spot"), (Patches, Cells, "Cell")):
            with self.subTest(layer=alias):
                sample = sample_class(bins.sample.gdf.copy(), bins.sample.adata.copy())
                args = (1., figure(), sample) if alias == "Spot" else (figure(), sample)
                layer = component_class(*args, alias=alias)
                layer.add_parameter(attr="alpha", colname="g0")
                layer.add_parameter(attr="color", static=True, val="#6baed6")
                self.assertIn("Deconvolution", [t.title for t in layer.get_widgets(Document()).tabs])

    def test_known_mixtures_dense_sparse_and_empty(self):
        for matrix in (np.array([[30., 10.], [0., 20.], [0., 0.]]),
                       sparse.csr_matrix([[30., 10.], [0., 20.], [0., 0.]])):
            target = ad.AnnData(matrix, var=pd.DataFrame(index=["g0", "g1"]))
            fractions, genes = deconvolve(target, reference(), "cell_type", batch_size=1)
            np.testing.assert_allclose(fractions.loc[:, ["A", "B"]].iloc[:2], [[.75, .25], [0, 1]])
            self.assertEqual(fractions.most_abundant.iloc[:2].tolist(), ["A", "B"])
            self.assertTrue(fractions.iloc[2].isna().all())
            self.assertEqual(genes, ["g0", "g1"])
            self.assertEqual(target.obsm, {})

    def test_invalid_inputs_fail_clearly(self):
        target = ad.AnnData(np.ones((2, 2)), var=pd.DataFrame(index=["g0", "g1"]))
        with self.assertRaisesRegex(ValueError, "cell-type column"):
            deconvolve(target, reference(), "missing")
        target.X[0, 0] = -1
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            deconvolve(target, reference(), "cell_type")
        target.X[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            deconvolve(target, reference(), "cell_type")
        target.var_names = ["other", "genes"]
        with self.assertRaisesRegex(ValueError, "shared genes"):
            deconvolve(target, reference(), "cell_type")
        ref = reference()
        ref.obs.loc["a", "cell_type"] = None
        with self.assertRaisesRegex(ValueError, "label"):
            deconvolve(reference(), ref, "cell_type")
        ref = reference()
        ref.X = sparse.csr_matrix(np.ones((3, 2)))
        with self.assertRaisesRegex(ValueError, "indistinguishable"):
            deconvolve(reference(), ref, "cell_type")

    def test_results_persist_and_bins_average_fractions(self):
        layer, _ = make_layer()
        result = deconvolve(layer.sample.adata, reference(), "cell_type")
        manager = layer.deconvolution_manager
        manager.apply_result(result, "cell_type")
        column = manager.columns[0]
        data = layer.sample.fetch_data(resolution_level="agg:2", colnames=[column, "g0"])
        self.assertTrue(data[column].between(0, 1).all())
        np.testing.assert_allclose(data[column].iloc[0], result[0].iloc[:2, 0].mean())
        self.assertEqual(data.g0.sum(), layer.sample.adata.X[:, 0].sum())
        self.assertIn(column, layer.autocomplete_input.completions)
        manager.apply_result(result, "cell_type")
        self.assertNotEqual(column, manager.columns[0])
        self.assertIn(column, layer.sample.adata.obs)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "results.h5ad"
            layer.sample.adata.write_h5ad(path)
            restored = ad.read_h5ad(path)
            pd.testing.assert_frame_equal(restored.obsm["deconvolution"], result[0])
            self.assertEqual(restored.uns["deconvolution"]["method"], "nnls")

    def test_fraction_display_remains_numeric_with_dominant_type_labels(self):
        layer, _ = make_layer()
        manager = layer.deconvolution_manager
        manager.apply_result(deconvolve(layer.sample.adata, reference(), "cell_type"), "cell_type")
        self.assertEqual(manager.columns, ["deconv_1:A", "deconv_1:B"])
        for name in manager.columns:
            self.assertTrue(pd.api.types.is_float_dtype(layer.sample.adata.obs[name]))
            self.assertTrue(layer.sample.is_continous(name))
        dominant = "deconv_1:most_abundant"
        self.assertFalse(layer.sample.is_continous(dominant))
        self.assertNotIn(dominant, layer.sample.adata.uns["mean_obs_columns"])
        self.assertNotIn(dominant, layer.sample.adata.uns["deconvolution"]["fraction_columns"])
        self.assertEqual(layer.sample.adata.uns["deconvolution"]["cell_types"], ["A", "B"])
        layer.show_colname(dominant)
        self.assertFalse(layer.manager.get_param("color").is_static)
        layer.show_colname(manager.columns[0])
        self.assertTrue(layer.manager.get_param("color").is_static)
        self.assertEqual(layer.manager.get_param("alpha").colname, manager.columns[0])
        self.assertFalse(layer.manager.get_param("alpha").avoid_scaling)
        alphas = np.asarray(layer.manager.source.data["alpha"], dtype=float)
        self.assertTrue(np.isfinite(alphas).all())
        self.assertGreater(len(np.unique(alphas)), 1)

    def test_tab_run_and_display_without_path_controls(self):
        layer, _ = make_layer()
        layer.deconvolution_manager.reference = reference()
        tabs = layer.get_widgets(ImmediateDoc())
        panel = next(t.child for t in tabs.tabs if t.title == "Deconvolution")
        widgets = panel.children
        self.assertFalse(any(isinstance(w, TextInput) for w in widgets))
        run, = [w for w in widgets if isinstance(w, Button)]
        labels = next(w for w in widgets if isinstance(w, Select) and w.title.startswith("Reference cell-type"))
        display = next(w for w in widgets if isinstance(w, Select) and w.title == "Display cell type")
        status = widgets[widgets.index(run) + 1]
        self.assertFalse(run.disabled)

        def wait_for(predicate):
            deadline = time.monotonic() + 10
            while not predicate() and time.monotonic() < deadline:
                time.sleep(.01)
            self.assertTrue(predicate(), status.text)

        self.assertEqual(labels.value, "cell_type")
        run._trigger_event(ButtonClick(run))
        wait_for(lambda: not run.disabled)
        self.assertIn("Deconvolved", status.text)
        self.assertFalse(display.disabled)
        self.assertEqual(layer.manager.get_param("alpha").colname, display.value)
        self.assertTrue(np.isfinite(layer.manager.source.data["alpha"]).all())


if __name__ == "__main__":
    unittest.main()
