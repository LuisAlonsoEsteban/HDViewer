"""Interaction cost: what is sent per bin, pruning, view masks and hidden layers."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from bokeh.models import AllIndices, BooleanFilter

from hd_viewer import settings
from hd_viewer.board.source.lazyloader import LazyLoader
from hd_viewer.board.source.parameter import CategoricalColor, StaticParameter
from hd_viewer.startup.sample.bins import block_mode
from tests.test_display_controls import make_layer


class BlockModeTests(unittest.TestCase):
    def test_majority_label_per_block_ignores_missing_and_breaks_ties_low(self):
        grid_x = np.array([0, 0, 0, 1, 1, 1, 1, 2])
        grid_y = np.zeros(8, dtype=int)
        values = np.array(["b", "a", "b", "c", "a", None, "c", None], dtype=object)
        modes = block_mode(grid_x, grid_y, values)
        self.assertEqual(modes.loc[(0, 0)], "b")
        self.assertEqual(modes.loc[(1, 0)], "c")
        self.assertTrue(pd.isna(modes.loc[(2, 0)]))           # only missing values
        tie = block_mode(np.zeros(2, int), np.zeros(2, int), np.array(["z", "a"], dtype=object))
        self.assertEqual(tie.iloc[0], "a")


class ColorCodeTests(unittest.TestCase):
    def test_categories_are_coded_against_a_growing_palette(self):
        param = CategoricalColor("kind", not_sequenced_val="#cccccc")
        param.category_colors = {"x": "#111111"}
        codes = param.get_values(pd.Series(["x", None, "y", "x"]))
        self.assertEqual(codes.dtype, np.uint16)
        self.assertEqual(param.decode(codes)[:2], ["#111111", "#cccccc"])
        before = list(codes)
        param.category_colors = {"x": "#222222"}               # recolour: palette only
        np.testing.assert_array_equal(param.get_values(pd.Series(["x", None, "y", "x"])), before)
        self.assertEqual(param.decode(codes)[0], "#222222")

    def test_static_colour_is_a_one_byte_code(self):
        param = StaticParameter("#abcdef", not_sequenced_val="#ffffff")
        codes = param.get_values([True, False])
        self.assertEqual(codes.dtype, np.uint8)
        self.assertEqual(param.decode(codes), ["#abcdef", "#ffffff"])


class StreamPayloadTests(unittest.TestCase):
    def test_bins_stream_compact_columns(self):
        layer = make_layer()
        data = layer.manager.source.data
        for name in ("width", "height", "angle"):
            self.assertNotIn(name, data)                       # glyph values instead
        for name in ("x", "y", "alpha"):
            self.assertEqual(np.asarray(data[name]).dtype, np.float32, name)
        self.assertEqual(np.asarray(data["color"]).dtype, np.uint8)
        self.assertEqual(layer.glyph_renderer.glyph.width, layer.sample.native_width)


class ViewFilterTests(unittest.TestCase):
    def test_no_mask_when_everything_is_visible_and_binary_mask_otherwise(self):
        layer = make_layer()
        layer._view = np.ones(layer.manager.n_geometries, dtype=bool)
        layer.change_glyph_view()
        self.assertIsInstance(layer.glyph_renderer.view.filter, AllIndices)
        mask = np.ones(layer.manager.n_geometries, dtype=bool)
        mask[:3] = False
        layer._view = mask
        layer.change_glyph_view()
        view_filter = layer.glyph_renderer.view.filter
        self.assertIsInstance(view_filter, BooleanFilter)
        self.assertIsInstance(view_filter.booleans, np.ndarray)
        mask = mask.copy()
        mask[5] = False
        layer._view = mask
        layer.change_glyph_view()
        self.assertIs(layer.glyph_renderer.view.filter, view_filter)   # mutated, not replaced
        np.testing.assert_array_equal(view_filter.booleans, mask)


class PruningTests(unittest.TestCase):
    def test_far_rows_are_kept_until_the_budget_is_exceeded(self):
        layer = make_layer()
        manager = layer.manager
        manager._last_fetch_rows = len(manager.indexes)
        self.assertFalse(manager.needs_pruning())
        manager._last_fetch_rows = len(manager.indexes) // (settings.PRUNE_ROW_BUDGET + 1)
        self.assertTrue(manager.needs_pruning())


class HiddenLayerTests(unittest.TestCase):
    def test_hidden_layers_drop_rows_and_skip_loads_until_shown(self):
        layer = make_layer()
        loader = LazyLoader(layer.manager, layer.figure)
        layer.loader = loader
        loaded = len(layer.manager.indexes)
        self.assertGreater(loaded, 0)

        layer.set_visibility(False)
        self.assertEqual(len(layer.manager.indexes), 0)
        self.assertEqual(len(layer.manager.source.data["x"]), 0)
        with patch.object(layer.sample, "fetch_data", wraps=layer.sample.fetch_data) as fetch:
            layer.manager.fetch_data()
            loader.force_reload()
            layer.show_colname("g1")                           # recolouring a hidden layer
            fetch.assert_not_called()

        layer.set_visibility(True)
        self.assertEqual(len(layer.manager.indexes), loaded)
        np.testing.assert_allclose(np.asarray(layer.manager.source.data["g1"], dtype=float),
                                   layer.sample.get_value(colname="g1").reindex(layer.manager.indexes).to_numpy(),
                                   rtol=1e-6)


if __name__ == "__main__":
    unittest.main()


class HoverNameTests(unittest.TestCase):
    def test_names_are_not_streamed_and_are_served_on_hover(self):
        import json
        from hd_viewer.board.components.layer.utils.hover_names import HoverNames
        layer = make_layer()
        manager, names = layer.manager, layer.hover_names
        self.assertNotIn("obs_name", manager.source.data)
        self.assertEqual(len(manager.loaded_obs_names), len(manager.source.data["x"]))
        self.assertEqual(layer.get_hover_tooltips(), HoverNames.TEMPLATE)
        self.assertIs(layer.hover.formatters["$index"], names.formatter)

        version = names.channel.version
        self.assertEqual(version, manager.row_version)
        names.channel.request = json.dumps(dict(seq=1, version=version, rows=[0, 2, 10_000]))
        response = json.loads(names.channel.response)
        self.assertEqual(response["rows"], [[0, manager.loaded_obs_names[0]], [2, manager.loaded_obs_names[2]]])

    def test_row_changes_invalidate_cached_names(self):
        import json
        layer = make_layer()
        manager, names = layer.manager, layer.hover_names
        before = names.channel.version
        manager.remove_indcs([manager.indexes[0]])            # positions shift
        self.assertGreater(names.channel.version, before)
        names.channel.request = json.dumps(dict(seq=2, version=before, rows=[0]))
        self.assertEqual(json.loads(names.channel.response)["rows"], [])   # stale request: no answer
        after = names.channel.version
        manager.clear_sources()
        self.assertGreater(names.channel.version, after)
