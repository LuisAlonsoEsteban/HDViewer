"""LayerRows without Bokeh: the operations a browser must apply for one layer."""
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.core.encoding import CategoricalColor, Scaler, StaticParameter
from hd_viewer.core.layer_state import LayerRows, RecordingSink
from hd_viewer.startup.sample import Bins


def make_rows():
    names = [f'b{i}' for i in range(64)]
    gdf = gpd.GeoDataFrame({'index': names}, index=names,
                           geometry=[box(i % 8, i // 8, i % 8 + 1, i // 8 + 1) for i in range(64)])
    adata = ad.AnnData(np.arange(128, dtype=float).reshape(64, 2),
                       obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['g0', 'g1']))
    adata.obs['kind'] = pd.Categorical(['a', 'b'] * 32)
    rows = LayerRows(Bins(gdf, adata, pixels_per_bin=1, bin_size=1))
    rows.add_parameter(Scaler(0., 1., 'g0', attr='alpha', not_sequenced_val=0., colname_hovered=True))
    rows.add_parameter(StaticParameter('#123456', not_sequenced_val='#888888'))
    rows.set_fetch_context(bbox=[-1, 9, -1, 9])
    rows.sink.take()
    return rows


class LayerRowsTests(unittest.TestCase):
    def test_a_load_appends_compact_columns_and_keeps_names_and_constants_server_side(self):
        rows = make_rows()
        rows.fetch_data()
        ops = rows.sink.take()
        self.assertEqual([op[0] for op in ops], ['append'])
        columns = ops[0][1]
        self.assertEqual(set(columns), {'x', 'y', 'alpha', 'color', 'g0'})
        self.assertEqual(columns['x'].dtype, np.float32)
        self.assertEqual(columns['color'].dtype, np.uint8)
        self.assertEqual(len(rows.loaded_obs_names), 64)
        self.assertEqual(rows.glyph_constants['width'], rows.sample.native_width)
        rows.fetch_data()                                     # same view: nothing new
        self.assertEqual(rows.sink.take(), [])

    def test_restyling_reuses_loaded_rows(self):
        rows = make_rows()
        rows.fetch_data()
        rows.sink.take()
        rows.update_param('color', new_val='#abcdef')         # palette only
        self.assertEqual(rows.sink.take(), [])
        self.assertEqual(set(rows.colors('color')), {'#abcdef'})
        rows.add_parameter(CategoricalColor('kind', attr='color', not_sequenced_val='#cccccc'))
        ops = rows.sink.take()
        self.assertTrue(ops and all(op[0] in ('set', 'replace') for op in ops))
        self.assertEqual(set(rows.colors('color')), set(rows.get_param('color').palette[1:]))
        version = rows.row_version
        rows.update_colname('g0', 'g1', attrs=('alpha',))
        self.assertEqual(rows.row_version, version)             # same rows, new values
        np.testing.assert_allclose(rows.sink.data['g1'], rows.sample.get_value(colname='g1').to_numpy())

    def test_pruning_and_hiding_change_row_positions(self):
        rows = make_rows()
        rows.fetch_data()
        rows.sink.take()
        version = rows.row_version
        rows.remove_further_data(0, 0, distance=3)
        ops = rows.sink.take()
        self.assertEqual(ops[-1][0], 'replace')
        self.assertLess(len(rows.loaded_obs_names), 64)
        self.assertEqual(len(rows.sink.data['x']), len(rows.loaded_obs_names))
        self.assertGreater(rows.row_version, version)
        rows.set_active(False)
        self.assertEqual(len(rows.indexes), 0)
        rows.fetch_data()
        self.assertEqual([op[0] for op in rows.sink.take()], ['replace'])     # only the clear


if __name__ == '__main__':
    unittest.main()
