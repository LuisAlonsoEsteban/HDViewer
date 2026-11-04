"""Numeric picker states and cluster display resolution constraints."""
import unittest
from unittest.mock import patch

from hd_viewer.startup.sample import Bins
import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import ColorPicker
from shapely.geometry import box

from hd_viewer.board.components.layer import Squares
from hd_viewer.board.source.lazyloader import LazyLoader


class ImmediateDoc:
    def add_next_tick_callback(self, callback):
        callback()


def make_layer():
    names = [f'b{i}' for i in range(64)]
    gdf = gpd.GeoDataFrame({'index': names}, index=names,
                          geometry=[box(i % 8, i // 8, i % 8 + 1, i // 8 + 1) for i in range(64)])
    adata = ad.AnnData(np.arange(128, dtype=float).reshape(64, 2),
                       obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['g0', 'g1']))
    adata.obs['score'] = np.linspace(0., 1., 64)
    adata.obs['kind'] = pd.Categorical(['a', 'b'] * 32)
    sample = Bins(gdf, adata, pixels_per_bin=1, bin_size=1)
    layer = Squares(figure(x_range=(-1, 9), y_range=(-1, 9)), sample, lazy_loading=True, alias='Bin')
    layer.add_parameter(attr='alpha', colname='g0', colname_hovered=True)
    layer.add_parameter(attr='color', static=True, val='#6baed6', not_sequenced_val='#b8b8b8')
    layer.manager.set_fetch_context(bbox=[-1, 9, -1, 9])
    layer.manager.fetch_data()
    layer.get_widgets(ImmediateDoc())
    return layer


class NumericPickerTests(unittest.TestCase):
    def test_gene_and_numeric_metadata_enable_picker(self):
        layer = make_layer()
        self.assertFalse(layer.colorpicker.disabled)
        for column in ('kind', 'g1', 'score', 'kind', 'g0'):
            layer.show_colname(column)
            self.assertEqual(layer.colorpicker.disabled, column == 'kind')

    def test_visibility_and_widget_rebuild_preserve_picker_state(self):
        layer = make_layer()
        layer.set_visibility(False)
        self.assertTrue(layer.colorpicker.disabled)
        layer.set_visibility(True)
        self.assertFalse(layer.colorpicker.disabled)
        layer.get_widgets(ImmediateDoc())
        self.assertFalse(layer.colorpicker.disabled)
        layer.show_colname('kind')
        layer.set_visibility(False)
        layer.set_visibility(True)
        self.assertTrue(layer.colorpicker.disabled)

    def test_numeric_picker_changes_color_without_changing_alpha(self):
        layer = make_layer()
        alpha = np.array(layer.manager.source.data['alpha'])
        layer.colorpicker.color = '#abcdef'
        self.assertEqual(set(layer.manager.source.data['color']), {'#abcdef'})
        np.testing.assert_array_equal(layer.manager.source.data['alpha'], alpha)
        layer.get_widgets(ImmediateDoc())
        self.assertEqual(layer.colorpicker.color, '#abcdef')

    def test_programmatic_color_change_cannot_mutate_categorical_display(self):
        layer = make_layer()
        layer.show_colname('kind')
        values = list(layer.manager.source.data['color'])
        layer.colorpicker.color = '#abcdef'
        self.assertEqual(list(layer.manager.source.data['color']), values)

    def test_extra_line_color_does_not_create_an_enabled_duplicate_picker(self):
        layer = make_layer()
        layer.add_parameter(attr='line_color', static=True, val='#abcdef', not_sequenced_val='#888888')
        layer.get_widgets(ImmediateDoc())
        pickers = [widget for widget in layer._all_widgets if isinstance(widget, ColorPicker)]
        self.assertEqual(pickers, [layer.colorpicker])
        self.assertFalse(layer.colorpicker.disabled)


class ClusterDisplayTests(unittest.TestCase):
    def clustered_layer(self):
        layer = make_layer()
        cm = layer.clustering_manager

        def fit(sample, *, key_added, **kwargs):
            sample.adata.obs[key_added] = pd.Categorical(['0', '1'] * (sample.adata.n_obs // 2))
            return key_added

        with patch.object(cm, '_cluster_on_sample', side_effect=fit):
            column = cm.run_clustering(resolution_level='agg:2', algorithm='kmeans')
        cm._build_cluster_state()
        cm.apply_colors_to_column(column)
        return layer, column

    def test_zoom_coarsens_cluster_display_but_not_below_analysis_level(self):
        layer, column = self.clustered_layer()
        labels = layer.sample.adata.obs[column].copy()
        analysis = layer.clustering_manager._get_active_analysis_sample()
        loader = LazyLoader(layer.manager, layer.figure)
        for span, expected in ((10, 2), (1000, 8), (10000, 12), (10, 2)):
            layer.figure.x_range.update(start=0, end=span)
            layer.figure.y_range.update(start=0, end=span)
            loader.force_reload()
            self.assertGreaterEqual(layer.sample.get_aggregation_factor(loader.get_resolution_level()), 2)
            self.assertEqual(layer.bin_info.metadata()['factor'], expected)
            self.assertEqual(layer.bin_info.mode_select.value, 'dynamic')
            self.assertTrue(all(name.endswith(f':{expected}') for name in layer.manager.source.data['obs_name']))
            pd.testing.assert_series_equal(layer.sample.adata.obs[column], labels)
            self.assertIs(layer.clustering_manager._get_active_analysis_sample(), analysis)

    def test_coarser_labels_colors_and_selection_match_native_members(self):
        layer, column = self.clustered_layer()
        layer.manager.set_aggregation_factor(4)
        LazyLoader(layer.manager, layer.figure).force_reload()
        names = layer.manager.source.data['obs_name']
        mapping = layer.sample.aggregation_names(4)
        for i, name in enumerate(names):
            members = mapping.index[mapping == name]
            expected = layer.sample.adata.obs.loc[members, column].value_counts().idxmax()
            self.assertEqual(layer.manager.source.data[column][i], expected)
            self.assertEqual(layer.manager.source.data['color'][i], layer.clustering_manager._color_map[expected])
        layer.set_selected_obs_names([layer.sample.adata.obs_names[0]])
        selected_name = mapping.iloc[0]
        self.assertEqual(layer.manager.source.selected.indices, [list(names).index(selected_name)])

    def test_switching_to_gene_restores_finer_levels_and_numeric_picker_state(self):
        layer, _ = self.clustered_layer()
        self.assertTrue(layer.colorpicker.disabled)
        self.assertNotIn('1', [value for value, _ in layer.bin_info.mode_select.options])
        layer.show_colname('g0')
        self.assertFalse(layer.colorpicker.disabled)
        self.assertIsNone(layer.manager.minimum_aggregation_factor)
        self.assertIn('1', [value for value, _ in layer.bin_info.mode_select.options])


if __name__ == '__main__':
    unittest.main()
