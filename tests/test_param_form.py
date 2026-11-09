"""Extra glyph parameters (Parameters tab) of a layer: core.parameters through LayerView.

The HTTP side (listing, validation, figures, frames) is in test_api_parameters_export."""
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.api.views import LayerView
from hd_viewer.core import parameters
from hd_viewer.core.encoding import CategoricalColor, Scaler, StaticParameter
from hd_viewer.startup.sample import Bins
from tests.fixtures import COLOUR, show


def make_layer(n=6):
    """A LayerView of n 16 µm bins in a row (unit squares [i, i+1] x [0, 1]) with a
    categorical `region` and a numeric `score` column, its rows loaded."""
    index = pd.Index([f"b{i}" for i in range(n)], name="barcode")
    gdf = gpd.GeoDataFrame({"index": list(index)}, geometry=[box(i, 0, i + 1, 1) for i in range(n)],
                           index=index)
    adata = ad.AnnData(np.arange(n * 2, dtype=float).reshape(n, 2),
                       obs=pd.DataFrame(index=index.astype(str)), var=pd.DataFrame(index=["g0", "g1"]))
    adata.obs["region"] = pd.Categorical(["a", "b"] * (n // 2))
    adata.obs["score"] = np.linspace(0, 1, n)
    view = LayerView("bins", "bins", Bins(gdf, adata, pixels_per_bin=1, bin_size=16), "Bins", COLOUR)
    show(view, -1, n + 1, -1, 2)
    return view


def data(view, attr):
    return np.asarray(view.rows.sink.data[attr])


class ParametersTests(unittest.TestCase):
    def test_only_relevant_columns_are_offered(self):
        view = make_layer()
        self.assertEqual(list(parameters.EXTRA_ATTRS), ["fill_alpha", "line_width", "line_color", "line_alpha"])
        widths = parameters.column_choices(view.sample, "line_width")
        self.assertIn("g0", widths)
        self.assertIn("score", widths)
        self.assertNotIn("region", widths)
        self.assertEqual(parameters.column_choices(view.sample, "line_color"), ["region"])
        self.assertEqual([entry["color"] for entry in view.describe_parameters()], [False, False, True, False])

    def test_add_edit_and_remove(self):
        view = make_layer()
        described = view.set_parameter("line_width", mode="static", value=2.5)
        self.assertEqual((described["mode"], described["value"]), ("static", 2.5))
        param = view.rows.get_param("line_width")
        self.assertIsInstance(param, StaticParameter)
        np.testing.assert_array_equal(data(view, "line_width"), [2.5] * 6)
        self.assertEqual(next(e["current"] for e in view.describe_parameters() if e["attr"] == "line_width"),
                         described)

        # Edit: saved as data-driven
        view.set_parameter("line_width", mode="data", column="score", minv=1.0, maxv=3.0)
        param = view.rows.get_param("line_width")
        self.assertIsInstance(param, Scaler)
        self.assertEqual((param.colname, param.minv, param.maxv), ("score", 1.0, 3.0))
        self.assertEqual(data(view, "line_width")[0], 1.0)
        self.assertEqual(data(view, "line_width")[-1], 3.0)
        current = next(e["current"] for e in view.describe_parameters() if e["attr"] == "line_width")
        self.assertEqual((current["mode"], current["column"], current["min"], current["max"]),
                         ("data", "score", 1.0, 3.0))

        view.remove_parameter("line_width")
        self.assertIsNone(view.rows.get_param("line_width"))
        self.assertNotIn("line_width", view.rows.sink.data)
        self.assertEqual(view.extra_parameters(), [])

    def test_categorical_line_colour_and_validation(self):
        view = make_layer()
        with self.assertRaisesRegex(ValueError, "categorical metadata"):
            view.set_parameter("line_color", mode="data")               # no column chosen
        with self.assertRaisesRegex(ValueError, "categorical metadata"):
            view.set_parameter("line_color", mode="data", column="score")
        self.assertIsNone(view.rows.get_param("line_color"))
        with self.assertRaisesRegex(ValueError, "between"):
            view.set_parameter("line_width", mode="static", value=1e6)
        with self.assertRaisesRegex(ValueError, "numeric"):
            view.set_parameter("line_width", mode="data", column="region")

        view.set_parameter("line_color", mode="data", column="region")
        self.assertIsInstance(view.rows.get_param("line_color"), CategoricalColor)
        self.assertEqual(len(set(data(view, "line_color"))), 2)

    def test_changing_the_colour_gene_keeps_data_parameters(self):
        view = make_layer()
        # Line width from the very gene the colour is showing (g0)
        self.assertEqual(view.rows.get_param("alpha").colname, "g0")
        view.set_parameter("line_width", mode="data", column="g0")
        before = list(data(view, "line_width"))

        view.colour_by_column("g1")

        self.assertEqual(view.rows.get_param("alpha").colname, "g1")
        self.assertEqual(view.rows.get_param("line_width").colname, "g0")
        self.assertEqual(list(data(view, "line_width")), before)

    def test_colour_and_alpha_are_not_offered(self):
        view = make_layer()
        self.assertFalse({"color", "alpha"} & set(parameters.EXTRA_ATTRS))
        self.assertEqual(view.extra_parameters(), [])
        for attr in ("color", "alpha"):
            with self.assertRaises(ValueError):
                view.set_parameter(attr, mode="static", value="#ffffff")
            with self.assertRaises(ValueError):
                view.remove_parameter(attr)


if __name__ == "__main__":
    unittest.main()
