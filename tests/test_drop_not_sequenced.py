"""Layers hold only geometries with counts: bins, pseudospots and cells without
expression are dropped when the sample is built (core.workspace.build_sample)."""
import unittest
from types import SimpleNamespace

import numpy as np

from hd_viewer.core.workspace import build_sample, sequenced_geometries
from tests.fixtures import expression, grid_frames, make_view, show


def b2s_stub():
    """The parts of a Bins2Spots that build_sample reads (no spatial transform)."""
    converter = SimpleNamespace(VisiumV2=SimpleNamespace(SPOT_DIAMETER=55.0),
                                visiumHD=SimpleNamespace(pixels_per_bin=1),
                                get_pixels=lambda um: um / 2)
    return SimpleNamespace(spatial_transform=None, converter=converter, bins_size=2,
                           tile_scale_xfactor=1.0, tile_scale_yfactor=1.0)


class DropNotSequencedTests(unittest.TestCase):
    def test_every_kind_keeps_only_geometries_with_counts(self):
        for kind, prefix in (("bins", "b"), ("spots", "s"), ("cells", "c")):
            with self.subTest(kind=kind):
                gdf, names = grid_frames(4, prefix=prefix, not_sequenced=3)     # 16 with counts + 3 without
                sequenced = names[::2]                                         # and half of the grid empty
                sample = build_sample(kind, gdf, expression(sequenced), b2s_stub())
                self.assertEqual(sorted(sample.gdf["index"]), sorted(sequenced))
                self.assertEqual(sample.n_geometries, sample.adata.n_obs)
                self.assertEqual(len(gdf), 19)                                  # the input is not modified

    def test_layers_load_only_sequenced_rows_also_when_aggregated(self):
        gdf, names = grid_frames(8, prefix="b", not_sequenced=2)
        sample = build_sample("bins", gdf, expression(names[: len(names) // 2]), b2s_stub())   # 4 of 8 rows
        for factor in (1, 2):
            with self.subTest(factor=factor):
                view = make_view(sample, kind="bins")
                view.set_aggregation(factor)
                show(view, -1, 9, -1, 9)
                self.assertTrue(len(view.rows.indexes))
                self.assertTrue(np.asarray(view.rows.is_sequenced, dtype=bool).all())
                self.assertIsNone(view.visibility())                            # nothing to hide

    def test_helper_is_a_no_op_when_everything_is_sequenced(self):
        gdf, names = grid_frames(3)
        self.assertIs(sequenced_geometries(gdf, expression(names)), gdf)
        self.assertIs(sequenced_geometries(gdf, None), gdf)


if __name__ == "__main__":
    unittest.main()
