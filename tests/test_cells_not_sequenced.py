"""Segmented cells without expression: kept in cells layers, hidden until
"Show not sequenced" is ticked, styled grey, named by their cell id on hover, and
left out of every analysis and of the UMIs/µm² area."""
import unittest

import numpy as np
import pandas as pd

from hd_viewer.core import summary
from hd_viewer.core.groups import GROUPS_COL
from tests.fixtures import column, loaded, make_view, show


def cells_view(not_sequenced=2):
    view = make_view(kind="cells", n=4, not_sequenced=not_sequenced)
    show(view, -1, 9, -1, 9)
    return view


def unsequenced(view):
    return pd.isna(np.asarray(loaded(view), dtype=object))


class CellsNotSequencedTests(unittest.TestCase):
    def test_loaded_but_hidden_by_default(self):
        view = cells_view()
        missing = unsequenced(view)
        self.assertEqual(missing.sum(), 2)                       # fetched with the others
        self.assertFalse(view.rows.show_not_sequenced)
        self.assertFalse(view.describe()["show_not_sequenced"])
        mask = view.visibility()
        self.assertIsNotNone(mask)
        np.testing.assert_array_equal(mask.astype(bool), ~missing)

    def test_toggle_shows_them_and_survives_reload_and_recolouring(self):
        view = cells_view()
        view.set_show_not_sequenced(True)
        self.assertIsNone(view.visibility())                     # every row drawn
        view.rows.clear_sources()
        view.rows.fetch_data()
        self.assertIsNone(view.visibility())
        view.set_show_not_sequenced(False)
        view.colour_by_column("g1")
        self.assertEqual((~view.visibility().astype(bool)).sum(), 2)
        view.rows.clear_sources()                                # reloaded rows keep the rule
        view.rows.fetch_data()
        np.testing.assert_array_equal(view.visibility().astype(bool), ~unsequenced(view))

    def test_styled_grey(self):
        view = cells_view()
        view.set_show_not_sequenced(True)
        alpha = np.asarray(column(view, "alpha"), dtype=float)
        missing = unsequenced(view)
        self.assertTrue(np.allclose(alpha[missing], 0.5))        # cells' not_sequenced_alpha
        self.assertEqual(set(np.asarray(column(view, "color"))[missing]), {0})   # not-sequenced colour code

    def test_groups_mode_follows_only_the_toggle(self):
        view = cells_view()
        obs_names = view.sample.adata.obs_names
        groups = view.groups
        groups.add_fixed_group("A", obs_names == obs_names[0])
        groups.refresh_active_labels()
        view.colour_by_categories(GROUPS_COL, groups._color_map, "groups")
        missing = unsequenced(view)
        self.assertFalse(view.visibility().astype(bool)[missing].any())
        view.set_show_not_sequenced(True)
        shown = view.visibility().astype(bool)
        self.assertTrue(shown[missing].all())
        self.assertEqual(shown[~missing].sum(), 1)               # only the group's member

    def test_hover_names_the_cell_and_selection_ignores_it(self):
        view = cells_view()
        missing = np.flatnonzero(unsequenced(view))
        names = dict(view.names(view.rows.row_version, [int(i) for i in missing]))
        self.assertTrue(all(name.startswith("not_sequenced") and "no expression" in name for name in names.values()))
        self.assertEqual(view.select(view.rows.row_version, [int(i) for i in missing]), 0)

    def test_summary_counts_and_area_exclude_them(self):
        view = cells_view(not_sequenced=3)
        stats = summary.general_stats(view.sample, "cells")
        self.assertEqual((stats["observations"], stats["segmented"], stats["not_sequenced"]), (16, 19, 3))
        self.assertEqual(stats["pct_sequenced"], round(100 * 16 / 19, 1))
        self.assertAlmostEqual(view.sample.get_unit_area_um2(), 1.0)   # unit squares only, not the 0.2 x 0.2 extras
        self.assertNotIn("segmented", summary.general_stats(view.sample, "bins"))

    def test_only_cells_layers_have_the_toggle(self):
        view = make_view(kind="bins")
        self.assertIsNone(view.describe()["show_not_sequenced"])
        with self.assertRaises(ValueError):
            view.set_show_not_sequenced(True)


if __name__ == "__main__":
    unittest.main()
