"""Interaction cost: what is sent per bin, pruning, view masks and hidden layers."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from hd_viewer import settings
from hd_viewer.api import frames
from hd_viewer.core.encoding import CategoricalColor, StaticParameter
from hd_viewer.startup.sample.bins import block_mode
from tests.fixtures import loaded, show
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
        view = make_layer()
        data = view.rows.sink.data
        for name in ("width", "height", "angle"):
            self.assertNotIn(name, data)                       # glyph values instead
        for name in ("x", "y", "alpha"):
            self.assertEqual(np.asarray(data[name]).dtype, np.float32, name)
        self.assertEqual(np.asarray(data["color"]).dtype, np.uint8)
        header = view.header([])
        self.assertEqual(view.rows.glyph_constants["width"], view.sample.native_width)
        self.assertEqual(header["constants"]["width"], view.sample.native_width)


class ViewFilterTests(unittest.TestCase):
    def test_no_mask_when_everything_is_visible_and_binary_mask_otherwise(self):
        view = make_layer()
        n = view.rows.n_geometries
        self.assertIsNone(view.visibility())
        header, _, extra = view.frame_parts([])
        self.assertFalse(header["filtered"])
        self.assertIsNone(extra)
        view.colour_by_column("kind")
        mask = np.ones(n, dtype=bool)
        mask[:3] = False
        view.rows.view_filter = mask
        visible = view.visibility()
        self.assertIsInstance(visible, np.ndarray)
        self.assertEqual(visible.dtype, np.uint8)               # one byte per row
        np.testing.assert_array_equal(visible, mask)
        header, ops, extra = view.frame_parts([])
        self.assertTrue(header["filtered"])
        np.testing.assert_array_equal(extra["visible"], mask)
        _, _, decoded = frames.decode(frames.encode(header, ops, extra))
        np.testing.assert_array_equal(np.asarray(decoded["visible"], dtype=bool), mask)


class PruningTests(unittest.TestCase):
    def test_far_rows_are_kept_until_the_budget_is_exceeded(self):
        rows = make_layer().rows
        rows._last_fetch_rows = len(rows.indexes)
        self.assertFalse(rows.needs_pruning())
        rows._last_fetch_rows = len(rows.indexes) // (settings.PRUNE_ROW_BUDGET + 1)
        self.assertTrue(rows.needs_pruning())


class HiddenLayerTests(unittest.TestCase):
    def test_hidden_layers_drop_rows_and_skip_loads_until_shown(self):
        view = make_layer()
        n_loaded = len(view.rows.indexes)
        self.assertGreater(n_loaded, 0)

        view.set_visible(False)
        self.assertEqual(len(view.rows.indexes), 0)
        self.assertEqual(len(view.rows.sink.data["x"]), 0)
        with patch.object(view.sample, "fetch_data", wraps=view.sample.fetch_data) as fetch:
            view.rows.fetch_data()
            show(view, -1, 9, -1, 9)
            view.colour_by_column("g1")                         # recolouring a hidden layer
            fetch.assert_not_called()

        view.set_visible(True)
        show(view, -1, 9, -1, 9)                                # the browser reloads its view
        self.assertEqual(len(view.rows.indexes), n_loaded)
        np.testing.assert_allclose(np.asarray(view.rows.sink.data["g1"], dtype=float),
                                   view.sample.get_value(colname="g1").reindex(view.rows.indexes).to_numpy(),
                                   rtol=1e-6)


class HoverNameTests(unittest.TestCase):
    def test_names_are_not_streamed_and_are_served_on_hover(self):
        view = make_layer()
        rows = view.rows
        self.assertNotIn("obs_name", rows.sink.data)
        self.assertEqual(len(rows.loaded_obs_names), len(rows.sink.data["x"]))
        self.assertEqual(view.tooltips()[0], ["Barcode", "$index{custom}"])
        version = rows.row_version
        self.assertEqual(view.header([])["row_version"], version)
        self.assertEqual(view.names(version, [0, 2, 10_000]),
                         [[0, rows.loaded_obs_names[0]], [2, rows.loaded_obs_names[2]]])

    def test_row_changes_invalidate_cached_names(self):
        view = make_layer()
        rows = view.rows
        before = rows.row_version
        rows.remove_indcs([rows.indexes[0]])                    # positions shift
        self.assertGreater(rows.row_version, before)
        self.assertEqual(view.names(before, [0]), [])           # stale request: no answer
        after = rows.row_version
        rows.clear_sources()
        self.assertGreater(rows.row_version, after)


if __name__ == "__main__":
    unittest.main()
