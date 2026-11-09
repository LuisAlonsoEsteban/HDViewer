"""Deconvolution of a layer: the single-cell reference (uploads, example paths,
annotations), NNLS, publishing results into the sample and the per-layer API."""
import os
import tempfile
import time
import unittest
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from fastapi.testclient import TestClient

from hd_viewer.analysis.deconvolution import deconvolve
from hd_viewer.api.sessions import COOKIE
from hd_viewer.api.views import build_views
from hd_viewer.core import deconvolution as core
from hd_viewer.core.pipeline import (LoadConfig, LoadedDataset, LoadPipeline,
                                     prepare_segmentation_image, read_deconvolution_reference)
from hd_viewer.core.session import session_upload_dir
from hd_viewer.startup import datasets
from hd_viewer.startup.data.singlecell import read_reference
from hd_viewer.startup.sample import Cells, Spots
from tests.test_param_form import make_layer
from tests.test_startup_uploads import UploadApiTestCase

# tests.test_core_pipeline (hence tests.test_api) imports tenx_reference_bytes from here:
# import neither of them at module level (tests.test_startup_uploads does not either).


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


def annotated_reference_file(folder, name="reference.h5ad"):
    """An h5ad reference with `cell_type` and `donor` annotations."""
    path = Path(folder) / name
    uploaded = reference()
    uploaded.obs["donor"] = ["d1", "d1", "d2"]
    uploaded.write_h5ad(path)
    return path


class DeconvolutionApiTestCase(UploadApiTestCase):
    """A session with the local runner (NNLS runs for real)."""

    def setUp(self):
        super().setUp()
        environment = patch.dict(os.environ, {"VHD_DECONV_RUNNER": "local"})
        environment.start()
        self.addCleanup(environment.stop)
        core.default_service.cache_clear()
        self.addCleanup(core.default_service.cache_clear)

    def install(self, samples, ref, base_dir, reference_name="reference.h5ad", label_column="", dataset=None):
        """Make `samples` ({id: (kind, sample, title)}) the session's loaded layers, as
        routers.load does once LoadPipeline returns."""
        session = self.session
        session.loaded = LoadedDataset(dataset, None, None, None, samples=dict(samples), reference=ref,
                                       reference_name=reference_name, reference_label_column=label_column,
                                       base_dir=Path(base_dir))
        session.views = build_views(session.loaded)
        return session

    def state(self, layer="bins"):
        response = self.client.get(f"/api/layers/{layer}/deconvolution")
        self.assertEqual(response.status_code, 200, response.text)
        return response.json()

    def run_method(self, layer="bins", method="nnls", label_column="cell_type", params=None, expect="done"):
        started = self.client.post(f"/api/layers/{layer}/deconvolution",
                                   json=dict(method=method, label_column=label_column, params=params or {}))
        self.assertEqual(started.status_code, 202, started.text)
        for _ in range(1200):
            state = self.state(layer)
            if state["job"]["state"] != "running":
                break
            time.sleep(.05)
        self.assertEqual(state["job"]["state"], expect, state["job"])
        return state


class ReferenceUploadTests(DeconvolutionApiTestCase):
    def test_uploaded_files_and_layer_results_use_each_sessions_folder(self):
        clients = [self.client, TestClient(self.client.app)]
        directories = []
        for client in clients:
            self.client = client
            self.client.get("/api/uploads")
            session = self.session
            self.assertEqual(session.directory, session_upload_dir(session.sid))
            self.assertEqual(session.directory.parent, Path(self.folder.name) / "uploads")
            self.upload("expression", "counts.h5", b"counts")
            saved = session.uploads["expression"].path
            self.assertEqual(saved, session.directory / "expression" / "counts.h5")
            self.assertEqual(saved.read_bytes(), b"counts")
            self.upload("positions", "positions.parquet")
            self.upload("image", "tissue.png")
            _, _, config = self.accepted()
            self.assertEqual(config.base_dir, session.directory)
            # LoadPipeline keeps base_dir; deconvolution results cache under it
            self.install({"bins": ("bins", make_layer().sample, "Bins")}, reference(), config.base_dir)
            self.run_method()
            log = Path(session.views["bins"].sample.adata.uns["deconvolution"]["log_path"])
            self.assertTrue(log.is_relative_to(session.directory / "deconvolution"))
            self.assertTrue(log.exists())
            directories.append(session.directory)
        self.assertNotEqual(*directories)

    def test_pending_reference_upload_blocks_load_until_it_arrives(self):
        self.upload_required()
        self.upload("reference_annotations", "labels.csv", b"barcode,cell_type\na,A\nb,A\nc,B\n")
        accepted, reply, _ = self.accepted()
        self.assertFalse(accepted)
        self.assertIn("not finished uploading", reply["detail"])
        self.assertEqual(self.status("reference")["status"], "empty")

        # The browser finishes the delayed reference upload chunk by chunk.
        self.upload("reference", "reference.h5", tenx_reference_bytes(), chunk_bytes=1000)
        self.assertEqual(self.status("reference")["status"], "complete")
        accepted, _, config = self.accepted()
        self.assertTrue(accepted)
        loaded = read_deconvolution_reference(config.reference_path, config.reference_annotations_path)
        self.assertEqual(loaded.obs.cell_type.tolist(), ["A", "A", "B"])

    def test_selected_reference_without_contents_is_not_silently_skipped(self):
        self.upload_required()
        self.assertEqual(self.start("reference", "reference.h5", 100).status_code, 200)
        accepted, reply, _ = self.accepted()
        self.assertFalse(accepted)
        self.assertIn("not finished uploading", reply["detail"])

    def test_tenx_h5_annotations_populate_columns(self):
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
                self.assertTrue(self.upload("reference", "counts.h5", tenx_reference_bytes(), 1000).json()["complete"])
                self.assertTrue(self.upload("reference_annotations", f"annotations.{format_name}", contents,
                                            1000).json()["complete"])
                fields = self.session.uploads
                loaded = read_deconvolution_reference(fields["reference"].path, fields["reference_annotations"].path)
                self.assertEqual(loaded.obs.cell_type.tolist(), ["A", "A", "B"])
                np.testing.assert_allclose(loaded.X.toarray(), reference().X.toarray())
                self.install({"bins": ("bins", make_layer().sample, "Bins")}, loaded, self.session.directory)
                state = self.state()
                self.assertEqual(state["unavailable"], "")
                self.assertEqual(state["columns"], ["cell_type", "donor"])
                self.assertEqual(state["label_column"], "cell_type")
                self.assertEqual(self.client.post("/api/layers/bins/deconvolution/check",
                                                  json=dict(method="nnls", label_column="cell_type")).status_code, 200)

    def test_counts_only_h5_explains_missing_annotations(self):
        self.install({"bins": ("bins", make_layer().sample, "Bins")}, read_reference(tenx_reference_bytes()),
                     self.session.directory)
        state = self.state()
        self.assertIn("no .obs annotation columns", state["unavailable"])
        self.assertEqual(state["columns"], [])
        body = dict(method="nnls", label_column="cell_type")
        self.assertEqual(self.client.post("/api/layers/bins/deconvolution/check", json=body).status_code, 422)
        self.assertEqual(self.client.post("/api/layers/bins/deconvolution", json=body).status_code, 422)
        self.assertEqual(self.state()["job"]["state"], "idle")
        with self.assertRaisesRegex(ValueError, "missing 2 reference cell barcodes"):
            read_reference(tenx_reference_bytes(), b"barcode,cell_type\na,A\n")
        with self.assertRaisesRegex(ValueError, "unique cell barcodes"):
            read_reference(tenx_reference_bytes(), b"barcode,cell_type\na,A\na,B\n")

    def test_available_on_spots_bins_and_cell_methods(self):
        bins = make_layer().sample
        samples = {"bins": ("bins", bins, "Bins"),
                   "spots": ("spots", Spots(bins.gdf.copy(), bins.adata.copy()), "Spots"),
                   "cells_a": ("cells", Cells(bins.gdf.copy(), bins.adata.copy()), "Cells · a")}
        self.install(samples, reference(), self.session.directory)
        for layer in samples:
            with self.subTest(layer=layer):
                state = self.state(layer)
                self.assertEqual(state["unavailable"], "")
                self.assertIn("nnls", [method["name"] for method in state["methods"]])
        self.assertEqual(self.client.get("/api/layers/nope/deconvolution").status_code, 404)


class ExampleReferenceLoadTests(DeconvolutionApiTestCase):
    def test_reference_is_read_once_and_reused_by_every_layer(self):
        data = Path(self.folder.name) / "data"
        path = annotated_reference_file(Path(self.folder.name))
        datasets.AVAILABLE_DATASETS["Tiny"] = dict(datasets.AVAILABLE_DATASETS["Tiny"], reference_path=str(path))
        with patch("anndata.read_h5ad", wraps=ad.read_h5ad) as read:
            self.load_example()
            session = self.session
            shared = session.loaded.reference
            self.assertEqual(shared.n_obs, 3)
            self.assertEqual(session.loaded.reference_name, "reference.h5ad")
            # A later layer of the same dataset uses the same reference
            session.loaded.samples["other"] = ("bins", make_layer().sample, "Other")
            session.views = build_views(session.loaded)
            for layer in ("bins", "other"):
                with self.subTest(layer=layer):
                    state = self.state(layer)
                    self.assertEqual(state["unavailable"], "")
                    self.assertEqual(state["reference"], dict(name="reference.h5ad", cells=3, genes=2))
                    self.assertEqual(state["label_column"], "cell_type")
                    self.assertEqual(state["columns"], ["cell_type", "donor"])
                    self.run_method(layer)
                    uns = session.views[layer].sample.adata.uns["deconvolution"]
                    self.assertEqual(uns["reference_name"], "reference.h5ad")
                    log = Path(uns["log_path"])
                    # Results cache next to the example dataset (its base folder)
                    self.assertTrue(log.is_relative_to(data / "deconvolution"))
                    self.assertTrue(log.exists())
            self.assertIs(session.loaded.reference, shared)
        self.assertEqual(read.call_count, 1)


class LoadPipelineReferenceTests(unittest.TestCase):
    def test_configured_dataset_uses_same_base_path_for_tiles_segmentation_and_deconvolution(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder) / "dataset"
            base.mkdir()
            image = base / "tissue.png"
            image.write_bytes(b"image")
            config = LoadConfig.from_example(dict(
                image_path=str(image), expression_path=str(base / "binned_outputs" / "square_008um" / "counts.h5"),
                positions_path=str(base / "binned_outputs" / "square_008um" / "positions.parquet"), bins_size=8))
            self.assertEqual(config.base_dir, base)
            recorded = []
            service = SimpleNamespace(prepare_image=lambda img, key, work_dir, bgr: recorded.append(work_dir)
                                      or work_dir / "image.tif")
            dataset = SimpleNamespace(img_path=str(image), img=np.zeros((2, 2, 3)), img_format="png",
                                      bins_size=8, pixels_per_bin=4, image_scale=1.0)
            prepared = prepare_segmentation_image(service, dataset, config.base_dir)
            self.assertEqual(recorded, [base / "segmentation"])
            self.assertEqual(prepared.work_dir, base / "segmentation")
            self.assertEqual(prepared.um_per_px, 2.0)

    def test_optional_or_invalid_reference(self):
        from tests.test_core_pipeline import write_dataset
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            example = write_dataset(root)
            loaded = LoadPipeline(LoadConfig.from_example(example, pseudospots=False)).run()
            self.assertIsNone(loaded.reference)
            self.assertIsNone(loaded.reference_name)
            self.assertEqual(loaded.base_dir, root)

            bad = root / "bad.h5ad"
            bad.write_bytes(b"not an h5ad file")
            config = LoadConfig.from_example(dict(example, reference_path=str(bad)), pseudospots=False)
            with patch("hd_viewer.core.pipeline.create_tiles") as tiles:
                with self.assertRaisesRegex(ValueError, "Cannot load single-cell reference bad.h5ad"):
                    LoadPipeline(config).run()
                tiles.assert_not_called()

    def test_parquet_annotations_validate_barcodes(self):
        for barcodes, message in ((["a"], "missing 2 reference cell barcodes"),
                                   (["a", "a"], "unique cell barcodes")):
            with self.subTest(barcodes=barcodes):
                buffer = BytesIO()
                pd.DataFrame({"barcode": barcodes, "cell_type": ["A"] * len(barcodes)}).to_parquet(buffer)
                with self.assertRaisesRegex(ValueError, message):
                    read_reference(tenx_reference_bytes(), buffer.getvalue())


class NNLSTests(unittest.TestCase):
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


class PublishedResultTests(unittest.TestCase):
    def test_results_persist_and_bins_average_fractions(self):
        view = make_layer()
        sample = view.sample
        result = deconvolve(sample.adata, reference(), "cell_type")
        self.assertEqual(core.publish_result(sample, result, "cell_type"), "1")
        column = core.runs(sample)["1"]["fraction_columns"][0]
        data = sample.fetch_data(resolution_level="agg:2", colnames=[column, "g0"])
        self.assertTrue(data[column].between(0, 1).all())
        np.testing.assert_allclose(data[column].iloc[0], result[0].iloc[:2, 0].mean())
        self.assertEqual(data.g0.sum(), sample.adata.X[:, 0].sum())
        self.assertIn(column, view.column_names("deconv"))
        self.assertEqual(core.publish_result(sample, result, "cell_type"), "2")
        self.assertNotEqual(column, core.runs(sample)["2"]["fraction_columns"][0])
        self.assertIn(column, sample.adata.obs)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "results.h5ad"
            sample.adata.write_h5ad(path)
            restored = ad.read_h5ad(path)
            pd.testing.assert_frame_equal(restored.obsm["deconvolution"], result[0])
            self.assertEqual(restored.uns["deconvolution"]["method"], "nnls")

    def test_fraction_display_remains_numeric_with_dominant_type_labels(self):
        view = make_layer()
        sample = view.sample
        core.publish_result(sample, deconvolve(sample.adata, reference(), "cell_type"), "cell_type")
        columns = core.runs(sample)["1"]["fraction_columns"]
        self.assertEqual(columns, ["deconv_1:A", "deconv_1:B"])
        for name in columns:
            self.assertTrue(pd.api.types.is_float_dtype(sample.adata.obs[name]))
            self.assertTrue(sample.is_continous(name))
        dominant = "deconv_1:most_abundant"
        self.assertFalse(sample.is_continous(dominant))
        self.assertNotIn(dominant, sample.adata.uns["mean_obs_columns"])
        self.assertNotIn(dominant, sample.adata.uns["deconvolution"]["fraction_columns"])
        self.assertEqual(sample.adata.uns["deconvolution"]["cell_types"], ["A", "B"])
        view.colour_by_column(dominant)
        self.assertFalse(view.rows.get_param("color").is_static)
        self.assertEqual(view.colour_by, dict(kind="category", name=dominant))
        view.colour_by_column(columns[0])
        self.assertTrue(view.rows.get_param("color").is_static)
        self.assertEqual(view.rows.get_param("alpha").colname, columns[0])
        self.assertFalse(view.rows.get_param("alpha").avoid_scaling)
        alphas = np.asarray(view.rows.sink.data["alpha"], dtype=float)
        self.assertTrue(np.isfinite(alphas).all())
        self.assertGreater(len(np.unique(alphas)), 1)


class RunAndDisplayTests(DeconvolutionApiTestCase):
    def test_run_then_colour_by_a_fraction(self):
        session = self.install({"bins": ("bins", make_layer().sample, "Bins")}, reference(), self.session.directory)
        state = self.state()
        self.assertEqual(state["runs"], [])
        state = self.run_method()
        self.assertIn("Deconvolved 6 observations into 2 cell types", state["job"]["message"])
        run, = state["runs"]
        self.assertEqual(run["fraction_columns"], ["deconv_1:A", "deconv_1:B"])
        self.client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=7, y0=-1, y1=2))
        response = self.client.patch("/api/layers/bins/style", json=dict(colour_by="deconv_1:A"))
        self.assertEqual(response.status_code, 200, response.text)
        view = session.views["bins"]
        self.assertEqual(view.rows.get_param("alpha").colname, "deconv_1:A")
        self.assertTrue(np.isfinite(np.asarray(view.rows.sink.data["alpha"], dtype=float)).all())


if __name__ == "__main__":
    unittest.main()
