import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.models import Div
from bokeh.events import ButtonClick
from bokeh.plotting import figure
from shapely.geometry import box

from hd_viewer.board.components.layer import Squares
from hd_viewer.board.source.parameter import CategoricalColor, Scaler, StaticParameter
from hd_viewer.startup.sample import Bins


def click(button):
    button._trigger_event(ButtonClick(button))


def make_layer(n=6):
    index = pd.Index([f"b{i}" for i in range(n)], name="barcode")
    gdf = gpd.GeoDataFrame({"index": list(index)}, geometry=[box(i, 0, i + 1, 1) for i in range(n)],
                           index=index)
    adata = ad.AnnData(np.arange(n * 2, dtype=float).reshape(n, 2),
                       obs=pd.DataFrame(index=index.astype(str)), var=pd.DataFrame(index=["g0", "g1"]))
    adata.obs["region"] = pd.Categorical(["a", "b"] * (n // 2))
    adata.obs["score"] = np.linspace(0, 1, n)
    layer = Squares(figure(), Bins(gdf, adata, pixels_per_bin=1, bin_size=16),
                    lazy_loading=True, alias="Bin")
    # As the app sets layers up
    layer.add_parameter(attr="alpha", colname="g0", colname_hovered=True)
    layer.add_parameter(attr="color", static=True, val="#6baed6", not_sequenced_val="#b8b8b8")
    layer.manager.bbox = (-1, n + 1, -1, 2)       # a loaded viewport
    layer.manager.fetch_data()
    layer.get_widgets(Document())
    return layer, layer.add_param_form


class ParamFormTests(unittest.TestCase):
    def test_only_the_relevant_inputs_are_shown(self):
        _, form = make_layer()
        self.assertEqual(form.attr_select.options, ["fill_alpha", "line_width", "line_color", "line_alpha"])
        form.attr_select.value = "line_width"
        visible = lambda: {name for name in ("value_color", "value_number", "colname_input", "range_row",
                                             "hover_checkbox", "ns_color", "ns_number")
                           if getattr(form, name).visible}
        self.assertEqual(visible(), {"value_number", "ns_number"})
        form.mode.active = 1                                   # from data
        self.assertEqual(visible(), {"colname_input", "range_row", "hover_checkbox", "ns_number"})
        self.assertIn("g0", form.colname_input.completions)
        self.assertIn("score", form.colname_input.completions)
        self.assertNotIn("region", form.colname_input.completions)

        form.attr_select.value = "line_color"                  # colours: categorical columns only
        self.assertEqual(visible(), {"colname_input", "hover_checkbox", "ns_color"})
        self.assertEqual(form.colname_input.completions, ["region"])
        form.mode.active = 0
        self.assertEqual(visible(), {"value_color", "ns_color"})

    def test_add_edit_and_remove(self):
        layer, form = make_layer()
        form.attr_select.value = "line_width"
        form.value_number.value = 2.5
        click(form.submit_button)
        param = layer.manager.get_param("line_width")
        self.assertIsInstance(param, StaticParameter)
        self.assertEqual(param.value, 2.5)
        # Drawn straight away: the viewport was reloaded with the new column
        self.assertEqual(layer.manager.source.data["line_width"], [2.5] * 6)
        self.assertEqual(layer.glyph_renderer.glyph.line_width, "line_width")   # drawn from the column
        self.assertNotIn("line_width", form.attr_select.options)
        self.assertEqual(len(form.list_column.children), 1)

        # Edit: the form is filled with the parameter, then saved as data-driven
        form._start_editing("line_width")
        self.assertTrue(form.attr_select.disabled)
        self.assertEqual(form.value_number.value, 2.5)
        form.mode.active = 1
        form.colname_input.value = "score"
        form.min_input.value, form.max_input.value = 1.0, 3.0
        click(form.submit_button)
        param = layer.manager.get_param("line_width")
        self.assertIsInstance(param, Scaler)
        self.assertEqual((param.colname, param.minv, param.maxv), ("score", 1.0, 3.0))
        self.assertEqual(layer.manager.source.data["line_width"][0], 1.0)
        self.assertEqual(layer.manager.source.data["line_width"][-1], 3.0)
        self.assertIn("updated", form.status.text)

        form._start_editing("line_width")                     # editing shows the saved values
        self.assertEqual((form.mode.active, form.colname_input.value), (1, "score"))

        click(form.list_column.children[0].children[2])       # ✕ remove
        self.assertIsNone(layer.manager.get_param("line_width"))
        self.assertNotIn("line_width", layer.manager.source.data)
        self.assertIn("line_width", form.attr_select.options)
        self.assertIn("No extra parameters", form.list_column.children[0].text)

    def test_categorical_line_colour_and_validation(self):
        layer, form = make_layer()
        form.attr_select.value = "line_color"
        form.mode.active = 1
        click(form.submit_button)                              # no column chosen
        self.assertIn("categorical metadata", form.status.text)
        self.assertIsNone(layer.manager.get_param("line_color"))

        form.colname_input.value = "region"
        click(form.submit_button)
        param = layer.manager.get_param("line_color")
        self.assertIsInstance(param, CategoricalColor)
        self.assertEqual(len(set(layer.manager.source.data["line_color"])), 2)

    def test_changing_the_colour_gene_keeps_data_parameters(self):
        layer, form = make_layer()
        # Line width from the very gene the Color tab is showing (g0)
        form.attr_select.value = "line_width"
        form.mode.active = 1
        form.colname_input.value = "g0"
        click(form.submit_button)
        before = list(layer.manager.source.data["line_width"])

        layer.show_colname("g1")                      # the Color tab's gene selector

        self.assertEqual(layer.manager.get_param("alpha").colname, "g1")
        self.assertEqual(layer.manager.get_param("line_width").colname, "g0")
        self.assertEqual(list(layer.manager.source.data["line_width"]), before)

    def test_colour_and_alpha_are_not_offered(self):
        _, form = make_layer()
        self.assertFalse({"color", "alpha"} & set(form.attr_select.options))
        self.assertEqual(form.extra_parameters(), [])
    def test_de_tab_uses_de_manager_widgets(self):
        layer, _ = make_layer()
        group_widgets = [Div(text="groups widget")]
        de_widgets = [Div(text="de widget")]

        layer.groups_manager.get_widgets = lambda doc: group_widgets
        layer.de_manager.get_widgets = lambda doc: de_widgets
        layer.clustering_manager.get_widgets = lambda doc, on_run=None, on_markers=None: []
        layer.deconvolution_manager.get_widgets = lambda doc: []
        layer.roi_manager.get_widgets = lambda doc: []

        layer.get_widgets(Document())

        de_tab = next(tab for tab in layer.tabs_widget.tabs if tab.title == "DE")
        self.assertIn(de_widgets[0], de_tab.child.children)
        self.assertNotIn(group_widgets[0], de_tab.child.children)

if __name__ == "__main__":
    unittest.main()
