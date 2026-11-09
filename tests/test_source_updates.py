"""Row updates preserve geometry, alignment, and compact browser operations."""
import threading
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

from hd_viewer.core.encoding import CategoricalColor, Scaler, StaticParameter
from hd_viewer.core.layer_state import LayerRows, source_columns


class Sample:
    def __init__(self, polygons=False):
        self.fetch_lock = threading.RLock()
        self.requests = []
        self.data = pd.DataFrame(dict(x=[1., 2., 3., 4.], y=[2., 3., 4., 5.],
                                     width=2., height=2., obs_name=['a', 'b', 'c', 'd'],
                                     gene=[0., np.nan, 5., 10.], other=[10., 20., 30., 40.]),
                                 index=['a', 'b', 'c', 'd'])
        if polygons:
            self.data['x'] = [[0., 1., 0.], [1., 2., 1.], [2., 3., 2.], [3., 4., 3.]]
            self.data['y'] = [[0., 0., 1.]] * 4
        self.adata = SimpleNamespace(obs_names=pd.Index(['a', 'c', 'd']))
        self.forced = False

    def fetch_data(self, colnames=None, bbox=None, return_coords=True, **kwargs):
        self.requests.append(dict(colnames=colnames, bbox=bbox, return_coords=return_coords, **kwargs))
        data = self.data if bbox is None else self.data.iloc[int(bbox[0]):int(bbox[1])]
        columns = ['obs_name', *(colnames or [])]
        if return_coords:
            columns = ['x', 'y', 'width', 'height', *columns]
        return data.loc[:, columns].copy()

    def get_range(self, name):
        return 0., 40.

    def is_continous(self, name):
        return True

    def get_forced_aggregation_factor(self, names):
        return 2 if self.forced and 'other' in names else None

    def get_resolution_level_for_factor(self, factor):
        return f'agg:{factor}'


class Source:
    """`manager.source.data` of the old tests: the sink's columns."""

    def __init__(self, sink):
        self.sink = sink

    @property
    def data(self):
        return self.sink.data


def make_manager(polygons=False, hovered=False, rows_class=LayerRows):
    sample = Sample(polygons)
    manager = rows_class(sample)
    manager.source = Source(manager.sink)
    manager.add_parameter(Scaler(0., 1., 'gene', attr='alpha', not_sequenced_val=0., colname_hovered=hovered))
    manager.add_parameter(StaticParameter('#123456', not_sequenced_val='#888888'))
    manager.bbox = [0, 4, 0, 4]
    manager.fetch_data()
    sample.requests.clear()
    manager.sink.take()
    return manager, sample, manager.sink


def column_ops(sink):
    """Operations sent since the last call, with the columns each one touches."""
    return [(op[0], {op[1]} if op[0] == 'set' else set(op[1])) for op in sink.take()]


class SourceUpdateTests(unittest.TestCase):
    def test_numeric_columns_are_owned_arrays_and_polygon_rings_remain_lists(self):
        sample = Sample(polygons=True)
        data = source_columns(sample.data)
        self.assertIsInstance(data['width'], np.ndarray)
        data['width'][0] = 99.
        self.assertEqual(sample.data.width.iloc[0], 2.)
        self.assertIsInstance(data['x'], list)
        self.assertEqual(data['x'][0], [0., 1., 0.])

    def test_static_replacement_only_sends_style(self):
        manager, sample, doc = make_manager(polygons=True)
        coordinates = manager.source.data['x']
        manager.add_parameter(StaticParameter('#abcdef', not_sequenced_val='#ffffff'))
        self.assertEqual(sample.requests, [])
        self.assertIs(manager.source.data['x'], coordinates)
        self.assertEqual(manager.colors('color'), ['#abcdef', '#ffffff', '#abcdef', '#abcdef'])
        # Geometry is never resent; unchanged colour codes are not resent either
        # (the new colours travel in the parameter's palette)
        ops = column_ops(doc)
        self.assertLessEqual(len(ops), 1)
        for kind, cols in ops:
            self.assertEqual((kind, cols), ('set', {'color'}))

    def test_non_hovered_scaling_reuses_raw_values(self):
        manager, sample, doc = make_manager()
        self.assertNotIn('gene', manager.source.data)
        manager.update_param('alpha', new_maxv=.8)
        self.assertEqual(sample.requests, [])
        np.testing.assert_allclose(manager.source.data['alpha'], [0., 0., .1, .2])
        self.assertEqual(column_ops(doc), [('set', {'alpha'})])

    def test_raw_cache_tracks_streamed_and_pruned_rows(self):
        manager, sample, _ = make_manager(polygons=True)
        manager.remove_indcs(['b', 'd'])
        self.assertEqual(manager.indexes, ['a', 'c'])
        self.assertEqual(manager.is_sequenced, [True, True])
        self.assertEqual(manager.source.data['x'], [[0., 1., 0.], [2., 3., 2.]])
        manager.update_param('alpha', new_maxv=.8)
        np.testing.assert_allclose(manager.source.data['alpha'], [0., .1])
        manager.fetch_data()
        self.assertEqual(manager.indexes, ['a', 'c', 'b', 'd'])
        manager.update_param('alpha', new_maxv=.8)
        np.testing.assert_allclose(manager.source.data['alpha'], [0., .1, 0., .2])
        self.assertEqual(len(manager.is_sequenced), 4)

    def test_gene_switch_fetches_values_without_replacing_geometry(self):
        manager, sample, doc = make_manager()
        coordinates = manager.source.data['x']
        sample.data = sample.data.iloc[::-1]
        manager.update_colname('gene', 'other')
        self.assertIs(manager.source.data['x'], coordinates)
        self.assertEqual(len(sample.requests), 1)
        self.assertFalse(sample.requests[0]['return_coords'])
        np.testing.assert_allclose(manager.source.data['alpha'], [.25, .5, .75, 1.])
        self.assertNotIn('gene', manager._value_cache)
        ops = column_ops(doc)
        self.assertTrue(ops)
        self.assertTrue(all(kind == 'set' and 'x' not in cols for kind, cols in ops))

    def test_incomplete_coverage_falls_back_to_full_fetch(self):
        manager, sample, _ = make_manager()
        manager.bbox = [2, 4, 0, 4]
        manager.update_colname('gene', 'other')
        self.assertEqual(manager.indexes, ['c', 'd'])
        self.assertEqual([request['return_coords'] for request in sample.requests], [False, True])
        np.testing.assert_allclose(manager.source.data['alpha'], [.75, 1.])

    def test_resolution_change_uses_full_fetch(self):
        manager, sample, _ = make_manager()
        sample.forced = True
        manager.update_colname('gene', 'other')
        self.assertEqual(len(sample.requests), 1)
        self.assertTrue(sample.requests[0]['return_coords'])
        self.assertEqual(sample.requests[0]['resolution_level'], 'agg:2')

    def test_hovered_raw_column_survives_shared_non_hovered_parameter(self):
        manager, sample, _ = make_manager(hovered=True)
        manager.add_parameter(Scaler(0., 2., 'gene', attr='line_width', not_sequenced_val=0.))
        self.assertIn('gene', manager.source.data)
        self.assertNotIn('gene', manager._value_cache)
        manager.remove_parameter('alpha')
        self.assertNotIn('gene', manager.source.data)
        self.assertIn('gene', manager._value_cache)
        sample.requests.clear()
        manager.update_param('line_width', new_maxv=4.)
        self.assertEqual(sample.requests, [])
        np.testing.assert_allclose(manager.source.data['line_width'], [0., 0., .5, 1.])

    def test_removing_parameter_drops_unused_cache_without_fetch(self):
        manager, sample, _ = make_manager()
        coordinates = manager.source.data['x']
        manager.remove_parameter('alpha')
        self.assertIs(manager.source.data['x'], coordinates)
        self.assertNotIn('alpha', manager.source.data)
        self.assertEqual(manager._value_cache, {})
        self.assertEqual(sample.requests, [])

    def test_categorical_values_remain_prunable_when_hover_column_is_removed(self):
        manager, sample, _ = make_manager()
        sample.data['kind'] = ['a', None, 'b', 'a']
        manager.add_parameter(CategoricalColor('kind', attr='line_color', colname_hovered=True))
        manager.add_parameter(CategoricalColor('kind', attr='color'))
        manager.remove_parameter('line_color')
        self.assertIsInstance(manager._value_cache['kind'], np.ndarray)
        manager.remove_indcs(['b'])
        self.assertEqual(manager._value_cache['kind'].tolist(), ['a', 'b', 'a'])

    def test_pruning_resynchronizes_selection_and_visibility(self):
        from tests.fixtures import loaded, make_view, show
        view = make_view(n=4)
        show(view, -1, 5, -1, 5)
        names = list(loaded(view))
        view.select(view.rows.row_version, [names.index('b1_1'), names.index('b2_2')])
        version = view.rows.row_version
        view.rows.remove_indcs(['b0_0', 'b1_1'])
        self.assertGreater(view.rows.row_version, version)        # old row positions are stale
        names = list(loaded(view))
        self.assertEqual(view.selected_positions(), [names.index('b2_2')])
        self.assertEqual(len(view.rows.view_filter), len(names))
        self.assertIsNone(view.visibility())                       # every remaining row is drawn
        with self.assertRaises(ValueError):
            view.select(version, [0])                              # positions of pruned rows

    def test_static_only_layers(self):
        sample = Sample()
        manager = LayerRows(sample)
        manager.add_parameter(StaticParameter(.5, attr='alpha'))
        manager.bbox = [0, 4, 0, 4]
        manager.fetch_data()
        manager.update_param('alpha', new_val=.7)
        np.testing.assert_allclose(manager.sink.data['alpha'], [.7] * 4)

    def test_real_samples_switch_values_without_replacing_geometry(self):
        import anndata as ad
        import geopandas as gpd
        from shapely.geometry import box
        from hd_viewer.startup.sample import Bins, Cells

        names = [f'b{i}' for i in range(16)]
        gdf = gpd.GeoDataFrame({'index': names}, index=names,
                              geometry=[box(i % 4, i // 4, i % 4 + 1, i // 4 + 1) for i in range(16)])
        adata = ad.AnnData(np.arange(32, dtype=float).reshape(16, 2),
                           obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['gene', 'other']))
        for kind, level in (('bins', None), ('bins', 'agg:2'), ('cells', None), ('cells', 'lod:0')):
            with self.subTest(kind=kind, level=level):
                sample = (Bins(gdf.copy(), adata.copy(), pixels_per_bin=1, bin_size=1)
                          if kind == 'bins' else Cells(gdf.copy(), adata.copy()))
                manager = LayerRows(sample)
                manager.source = Source(manager.sink)
                manager.add_parameter(Scaler(0., 1., 'gene', attr='alpha', not_sequenced_val=0.))
                manager.set_fetch_context(bbox=[-1, 5, -1, 5], resolution_level=level)
                manager.fetch_data()
                coordinates = manager.source.data['x']
                indexes = manager.indexes.copy()
                manager.update_colname('gene', 'other')
                self.assertIs(manager.source.data['x'], coordinates)
                self.assertEqual(manager.indexes, indexes)
                expected = sample.fetch_data(bbox=manager.bbox, colnames=['other'], resolution_level=level)
                np.testing.assert_allclose(manager._value_cache['other'], expected['other'])

    def test_value_only_fetch_preserves_real_observation_index(self):
        import anndata as ad
        import geopandas as gpd
        from shapely.geometry import box
        from hd_viewer.startup.sample import Cells

        names = ['a', 'b']
        gdf = gpd.GeoDataFrame({'index': names}, index=names,
                              geometry=[box(0, 0, 1, 1), box(2, 0, 3, 1)])
        adata = ad.AnnData(np.array([[2.], [7.]]), obs=pd.DataFrame(index=names),
                           var=pd.DataFrame(index=['gene']))
        sample = Cells(gdf, adata)
        fetched = sample.fetch_data(colnames=['gene'], return_coords=False)
        self.assertEqual(fetched.index.tolist(), names)
        np.testing.assert_array_equal(fetched.gene, [2., 7.])
        self.assertNotIn('x', fetched)


if __name__ == '__main__':
    unittest.main()
