"""Genes can be displayed as counts, normalized or log-normalized counts."""
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.startup.sample import Bins
from tests.fixtures import column, loaded
from tests.test_roi_manager import Api, make_layer


def expected(counts, library, scale):
    if scale == "counts":
        return counts
    normalized = np.where(library == 0, 0.0, counts * 1e4 / np.where(library == 0, 1, library))
    return np.log1p(normalized) if scale == "log1p" else normalized


class SampleScaleTests(unittest.TestCase):
    def make_bins(self):
        ids = [str(i) for i in range(64)]
        gdf = gpd.GeoDataFrame({"index": ids + ["unsequenced"]},
                               geometry=[box(x * 2, y * 2, x * 2 + 2, y * 2 + 2)
                                         for y in range(8) for x in range(8)] + [box(0, 0, 1, 1)],
                               index=pd.Index(ids + ["unsequenced"], name="barcode"))
        rng = np.random.default_rng(0)
        X = rng.poisson(3, size=(64, 3)).astype(float)
        X[5] = 0                                    # an empty observation
        adata = ad.AnnData(X, obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=["g0", "g1", "g2"]))
        return Bins(gdf, adata, pixels_per_bin=2), X

    def test_native_values_and_ranges_follow_the_scale(self):
        sample, X = self.make_bins()
        library = X.sum(axis=1)
        for scale in ("counts", "normalized", "log1p"):
            with self.subTest(scale=scale):
                sample.expression_scale = scale
                values = sample.get_value(colname="g0").to_numpy()
                np.testing.assert_allclose(values[:64], expected(X[:, 0], library, scale))
                self.assertTrue(np.isnan(values[64]))          # not sequenced stays NaN
                self.assertEqual(values[5], 0)                 # 0/0 shown as 0
                low, high = sample.get_range("g0")
                self.assertAlmostEqual(high, expected(X[:, 0], library, scale).max())
                # Analyses read counts regardless of the display
                np.testing.assert_allclose(sample.get_value(colname="g0", raw=True).to_numpy()[:64], X[:, 0])

    def test_aggregates_normalize_by_their_own_total_counts(self):
        sample, X = self.make_bins()
        sample.expression_scale = "counts"
        counts = sample.fetch_data(resolution_level="agg:2", colnames=["g0", "g1", "g2"])
        for scale in ("normalized", "log1p"):
            with self.subTest(scale=scale):
                sample.expression_scale = scale
                data = sample.fetch_data(resolution_level="agg:2", colnames=["g0"])
                library = counts[["g0", "g1", "g2"]].sum(axis=1).loc[data.index].to_numpy()
                np.testing.assert_allclose(
                    data["g0"].to_numpy(), expected(counts.loc[data.index, "g0"].to_numpy(), library, scale))

    def test_unknown_scale_is_rejected(self):
        sample, _ = self.make_bins()
        with self.assertRaises(ValueError):
            sample.expression_scale = "zscore"


class LayerScaleTests(unittest.TestCase):
    def test_scale_change_redraws_every_layer_type(self):
        for kind in ("bins", "spots", "cells"):
            with self.subTest(kind=kind):
                api = Api(self, kind)
                view = api.view
                self.assertEqual(view.describe()["scale"], "counts")
                X = np.asarray(view.sample.adata.X)
                library = dict(zip(view.sample.adata.obs_names, X.sum(axis=1)))
                g0 = dict(zip(view.sample.adata.obs_names, X[:, 0]))
                self.assertEqual(view.colour_by, dict(kind="gene", name="g0"))

                for scale in ("log1p", "normalized", "counts"):
                    response = api.patch("style", json=dict(scale=scale))
                    self.assertEqual(response.status_code, 200, response.text)
                    self.assertEqual(view.describe()["scale"], scale)
                    names = loaded(view)
                    shown = np.asarray(column(view, "g0"), dtype=float)
                    sequenced = [name in g0 for name in names]
                    self.assertEqual(sum(not s for s in sequenced), 1)          # the not-sequenced shape
                    want = expected(np.array([g0[n] for n in names if n in g0]),
                                    np.array([library[n] for n in names if n in g0]), scale)
                    np.testing.assert_allclose(shown[sequenced], want, rtol=1e-6)
                    alpha = view.rows.get_param("alpha")
                    self.assertAlmostEqual(alpha.max_counts, view.sample.get_range("g0")[1], places=5)
                self.assertEqual(api.patch("style", json=dict(scale="zscore")).status_code, 422)
                self.assertEqual(view.describe()["scale"], "counts")

    def test_tooltips_name_the_scale_of_the_coloured_gene(self):
        view = make_layer("bins")
        for scale, label in (("counts", "g0"), ("normalized", "g0 (normalized)"),
                             ("log1p", "g0 (log-normalized)")):
            view.set_scale(scale)
            self.assertIn([label, "@{g0}"], view.tooltips())

    def test_tooltips_name_the_scale_of_hovered_parameter_genes(self):
        view = make_layer("bins")
        view.set_parameter("line_width", mode="data", column="g1", hovered=True)
        for scale, label in (("counts", "g1"), ("normalized", "g1 (normalized)"),
                             ("log1p", "g1 (log-normalized)")):
            view.set_scale(scale)
            self.assertIn([label, "@{g1}"], view.tooltips())


if __name__ == "__main__":
    unittest.main()
