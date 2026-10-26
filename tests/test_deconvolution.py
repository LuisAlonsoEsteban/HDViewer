import tempfile
import time
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import Button, Select, TextInput
from bokeh.plotting import figure
from scipy import sparse

from src.deconvolution import deconvolve
from src.sample import Cells, Spots
from bokeh_server.components.glyphs import Circles, Patches
from tests.test_param_form import make_layer


def reference():
    # Deliberately reversed genes and unequal library sizes.
    return ad.AnnData(sparse.csr_matrix([[0., 100.], [0., 10.], [40., 0.]]),
                      obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["a", "b", "c"]),
                      var=pd.DataFrame(index=["g1", "g0"]))


class ImmediateDoc(Document):
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


class DeconvolutionTests(unittest.TestCase):
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
            np.testing.assert_allclose(fractions.iloc[:2], [[.75, .25], [0, 1]])
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

    def test_tab_load_run_and_display(self):
        layer, _ = make_layer()
        tabs = layer.get_widgets(ImmediateDoc())
        panel = next(t.child for t in tabs.tabs if t.title == "Deconvolution")
        widgets = panel.children
        path = next(w for w in widgets if isinstance(w, TextInput))
        load, run = [w for w in widgets if isinstance(w, Button)]
        labels, display = [w for w in widgets if isinstance(w, Select)]
        status = widgets[widgets.index(run) + 1]

        def wait_for(predicate):
            deadline = time.monotonic() + 10
            while not predicate() and time.monotonic() < deadline:
                time.sleep(.01)
            self.assertTrue(predicate(), status.text)

        with tempfile.TemporaryDirectory() as folder:
            ref_path = Path(folder) / "reference.h5ad"
            reference().write_h5ad(ref_path)
            path.value = str(ref_path)
            load._trigger_event(ButtonClick(load))
            wait_for(lambda: not load.disabled)
            self.assertEqual(labels.value, "cell_type")
            self.assertFalse(run.disabled)
            run._trigger_event(ButtonClick(run))
            wait_for(lambda: not run.disabled)
            self.assertIn("Deconvolved", status.text)
            self.assertFalse(display.disabled)
            self.assertEqual(layer.manager.get_param("alpha").colname, display.value)
            self.assertTrue(np.isfinite(layer.manager.source.data["alpha"]).all())


if __name__ == "__main__":
    unittest.main()
