"""Colour controls (continuous vs categorical) and cluster display resolution constraints."""
import unittest
from unittest.mock import patch

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.core.encoding import CategoricalColor, StaticParameter
from hd_viewer.startup.sample import Bins
from tests.fixtures import loaded, make_view, show


def make_bins():
    names = [f'b{i}' for i in range(64)]
    gdf = gpd.GeoDataFrame({'index': names}, index=names,
                           geometry=[box(i % 8, i // 8, i % 8 + 1, i // 8 + 1) for i in range(64)])
    adata = ad.AnnData(np.arange(128, dtype=float).reshape(64, 2),
                       obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['g0', 'g1']))
    adata.obs['score'] = np.linspace(0., 1., 64)
    adata.obs['kind'] = pd.Categorical(['a', 'b'] * 32)
    return Bins(gdf, adata, pixels_per_bin=1, bin_size=1)


def make_layer():
    """The view of an 8 x 8 bin layer coloured by g0, showing every bin."""
    view = make_view(make_bins(), colour='#6baed6')
    show(view, -1, 9, -1, 9)
    return view


def picker_enabled(view):
    """The colour picker applies to continuous values only (one colour, scaled opacity)."""
    return view.colour_by['kind'] in ('gene', 'metadata')


class NumericPickerTests(unittest.TestCase):
    def test_gene_and_numeric_metadata_enable_picker(self):
        view = make_layer()
        self.assertTrue(picker_enabled(view))
        for column in ('kind', 'g1', 'score', 'kind', 'g0'):
            view.colour_by_column(column)
            self.assertEqual(picker_enabled(view), column != 'kind')
            colour = view.rows.get_param('color')
            self.assertIsInstance(colour, CategoricalColor if column == 'kind' else StaticParameter)
            self.assertEqual(view.describe()['colour_by']['name'], column)

    def test_visibility_preserves_picker_state(self):
        view = make_layer()
        view.set_visible(False)
        view.set_visible(True)
        self.assertTrue(picker_enabled(view))
        view.colour_by_column('kind')
        view.set_visible(False)
        view.set_visible(True)
        show(view, -1, 9, -1, 9)
        self.assertFalse(picker_enabled(view))
        self.assertIsInstance(view.rows.get_param('color'), CategoricalColor)

    def test_numeric_picker_changes_color_without_changing_alpha(self):
        view = make_layer()
        alpha = np.array(view.rows.sink.data['alpha'])
        codes = np.array(view.rows.sink.data['color'])
        ops = view.set_colour('#abcdef')
        self.assertEqual(set(view.rows.colors('color')), {'#abcdef'})
        # Only the palette changes: the per-row colour codes are not resent
        self.assertEqual(ops, [])
        self.assertIn('#abcdef', view.header(ops)['palette'])
        np.testing.assert_array_equal(view.rows.sink.data['color'], codes)
        np.testing.assert_array_equal(view.rows.sink.data['alpha'], alpha)
        self.assertEqual(view.describe()['colour'], '#abcdef')

    def test_color_change_cannot_mutate_categorical_display(self):
        view = make_layer()
        view.colour_by_column('kind')
        values = list(view.rows.sink.data['color'])
        colours = view.rows.colors('color')
        view.set_colour('#abcdef')
        self.assertEqual(list(view.rows.sink.data['color']), values)
        self.assertEqual(view.rows.colors('color'), colours)
        self.assertNotIn('#abcdef', colours)
        view.colour_by_column('g0')                     # back to continuous: the new colour
        self.assertEqual(set(view.rows.colors('color')), {'#abcdef'})

    def test_extra_line_color_does_not_replace_the_layer_colour(self):
        view = make_layer()
        view.set_parameter('line_color', mode='static', value='#abcdef', not_sequenced='#888888')
        self.assertEqual(view.rows.get_param('color').value, '#6baed6')
        self.assertTrue(picker_enabled(view))
        self.assertEqual(set(view.rows.colors('color')), {'#6baed6'})
        self.assertEqual(set(view.rows.colors('line_color')), {'#abcdef'})


def clustered_layer():
    view = make_layer()
    clustering = view.clustering

    def fit(sample, *, key_added, **kwargs):
        sample.adata.obs[key_added] = pd.Categorical(['0', '1'] * (sample.adata.n_obs // 2))
        return key_added

    with patch.object(clustering, '_cluster_on_sample', side_effect=fit):
        column = clustering.compute(1.0, resolution_level='agg:2', algorithm='kmeans')
    clustering._build_cluster_state()
    view.colour_by_categories(column, clustering._color_map, 'clusters')
    return view, column


def zoom(view, span):
    return show(view, 0, span, 0, span)


class ClusterDisplayTests(unittest.TestCase):
    def test_zoom_coarsens_cluster_display_but_not_below_analysis_level(self):
        view, column = clustered_layer()
        labels = view.sample.adata.obs[column].copy()
        analysis = view.clustering._get_active_analysis_sample()
        for span, expected in ((10, 2), (1000, 8), (10000, 12), (10, 2)):
            ops = zoom(view, span)
            level = view.rows.effective_resolution_level
            self.assertGreaterEqual(view.sample.get_aggregation_factor(level), 2)
            self.assertEqual(view.sample.get_aggregation_factor(level), expected)
            self.assertEqual(view.header(ops)['level'], f'agg:{expected}')
            self.assertIsNone(view.describe()['fixed_aggregation'])     # still follows the zoom
            self.assertTrue(len(loaded(view)))
            self.assertTrue(all(name.endswith(f':{expected}') for name in loaded(view)))
            pd.testing.assert_series_equal(view.sample.adata.obs[column], labels)
            self.assertIs(view.clustering._get_active_analysis_sample(), analysis)

    def test_finer_fixed_level_than_the_analysis_is_refused(self):
        view, _ = clustered_layer()
        self.assertEqual(view.rows.minimum_aggregation_factor, 2)
        with self.assertRaises(ValueError):
            view.set_aggregation(1)

    def test_coarser_labels_colors_and_selection_match_native_members(self):
        view, column = clustered_layer()
        view.set_aggregation(4)
        zoom(view, 10)
        names = loaded(view)
        self.assertTrue(names)
        mapping = view.sample.aggregation_names(4)
        colours = view.rows.colors('color')
        for i, name in enumerate(names):
            members = mapping.index[mapping == name]
            expected = view.sample.adata.obs.loc[members, column].value_counts().idxmax()
            self.assertEqual(view.rows.sink.data[column][i], expected)
            self.assertEqual(colours[i], view.clustering._color_map[expected])
        view.selected = {view.sample.adata.obs_names[0]}
        selected_name = mapping.iloc[0]
        self.assertEqual(view.selected_positions(), [list(names).index(selected_name)])

    def test_switching_to_gene_restores_finer_levels_and_numeric_picker_state(self):
        view, _ = clustered_layer()
        self.assertFalse(picker_enabled(view))
        self.assertEqual(view.rows.minimum_aggregation_factor, 2)
        view.colour_by_column('g0')
        self.assertTrue(picker_enabled(view))
        self.assertIsNone(view.rows.minimum_aggregation_factor)
        view.set_aggregation(1)                               # native bins allowed again
        self.assertEqual(view.rows.effective_resolution_level, view.sample.get_resolution_level_for_factor(1))


if __name__ == '__main__':
    unittest.main()
