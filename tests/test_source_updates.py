"""Source updates preserve geometry, alignment, and compact Bokeh messages."""
import threading
import unittest
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.document.events import ColumnDataChangedEvent

from hd_viewer.board.source.manager import DynamicSourceManager, SourceManager, source_columns
from hd_viewer.board.source.parameter import CategoricalColor, Scaler, StaticParameter


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


def make_manager(polygons=False, hovered=False):
    sample = Sample(polygons)
    manager = DynamicSourceManager(sample)
    manager.add_parameter(Scaler(0., 1., 'gene', attr='alpha', not_sequenced_val=0., colname_hovered=hovered))
    manager.add_parameter(StaticParameter('#123456', not_sequenced_val='#888888'))
    manager.bbox = [0, 4, 0, 4]
    manager.fetch_data()
    sample.requests.clear()
    doc = Document()
    doc.add_root(manager.source)
    doc.hold('collect')
    return manager, sample, doc


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
        self.assertEqual(manager.source.data['color'], ['#abcdef', '#ffffff', '#abcdef', '#abcdef'])
        events = doc.callbacks._held_events
        self.assertEqual(len(events), 1)
        self.assertIsInstance(events[0], ColumnDataChangedEvent)
        self.assertEqual(set(events[0].cols), {'color'})

    def test_non_hovered_scaling_reuses_raw_values(self):
        manager, sample, doc = make_manager()
        self.assertNotIn('gene', manager.source.data)
        manager.update_param('alpha', new_maxv=.8)
        self.assertEqual(sample.requests, [])
        np.testing.assert_allclose(manager.source.data['alpha'], [0., 0., .1, .2])
        self.assertEqual(set(doc.callbacks._held_events[0].cols), {'alpha'})

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
        self.assertNotIn('x', doc.callbacks._held_events[0].cols)

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
        manager, _, _ = make_manager()
        manager.component = SimpleNamespace(show_not_sequenced=True,
                                            apply_selected_obs_names=Mock(),
                                            update_hover_tooltips=Mock(), change_glyph_view=Mock())
        manager.remove_indcs(['a', 'b'])
        manager.component.apply_selected_obs_names.assert_called_once()
        manager.component.change_glyph_view.assert_called_once()
        self.assertEqual(manager.component._view, [True, True])

    def test_hidden_unsequenced_mask_is_rebuilt_on_gene_switch(self):
        manager, _, _ = make_manager()
        manager.component = SimpleNamespace(show_not_sequenced=False,
                                            apply_selected_obs_names=Mock(),
                                            update_hover_tooltips=Mock(), change_glyph_view=Mock())
        manager.view_filter = np.array([True, False, True, True])
        manager.update_colname('gene', 'other')
        self.assertEqual(manager.component._view, [True] * 4)

    def test_static_only_layers_and_non_lazy_sources(self):
        sample = Sample()
        manager = DynamicSourceManager(sample)
        manager.add_parameter(StaticParameter(.5, attr='alpha'))
        manager.bbox = [0, 4, 0, 4]
        manager.fetch_data()
        manager.update_param('alpha', new_val=.7)
        np.testing.assert_allclose(manager.source.data['alpha'], [.7] * 4)
        static = SourceManager(sample)
        static.add_parameter(StaticParameter(.5, attr='alpha'))
        static.fetch_data()
        static.update_param('alpha', new_val=.7)
        np.testing.assert_allclose(static.source.data['alpha'], [.7] * 4)

    def test_real_samples_switch_values_without_replacing_geometry(self):
        # Initialize in the same order as the app (startup imports the UI).
        import hd_viewer.startup
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
                manager = DynamicSourceManager(sample)
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
        import hd_viewer.startup
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
