import threading
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import Div, HoverTool

from hd_viewer.board.operations.layer.metadata_comparison import compare_metadata, compare_metadata_to_rest, rest_mask, numeric_columns, build_metadata_plot
from hd_viewer.board.operations.layer.deconvolution import DeconvolutionResult
# Follow application import order before loading layer UI fixtures.
import hd_viewer.startup
from tests.test_roi_manager import make_layer, ImmediateDoc, wait_for


class MetadataSummaryTests(unittest.TestCase):
    def test_group_means_preserve_proportions_and_count_finite_values(self):
        frame = pd.DataFrame({'T': [.1, .3, .8, np.nan], 'B': [.9, .7, .2, np.inf], 'category': ['a'] * 4})
        masks = {'A': np.array([1, 1, 0, 0], dtype=bool), 'B': np.array([0, 0, 1, 1], dtype=bool)}
        before = frame.copy()
        result = compare_metadata(frame, masks, ['T', 'B'])
        values = result.summary.set_index(['group', 'column'])
        self.assertAlmostEqual(values.loc[('A', 'T'), 'value'], .2)
        self.assertAlmostEqual(values.loc[('B', 'T'), 'value'], .8)
        self.assertEqual(values.loc[('B', 'T'), 'n'], 1)
        self.assertEqual(values.loc[('B', 'T'), 'total'], 2)
        # Selecting one cell type keeps its fraction rather than scaling it to one.
        only = compare_metadata(frame, masks, ['T'])
        self.assertAlmostEqual(only.summary.value.iloc[0], .2)
        pd.testing.assert_frame_equal(frame, before)

    def test_overlapping_observations_and_median(self):
        frame = pd.DataFrame({'score': [0., 2., 10.]})
        masks = {'A': [True, True, False], 'B': [False, True, True]}
        result = compare_metadata(frame, masks, ['score'], 'median')
        self.assertEqual(result.summary.value.tolist(), [1., 6.])
        self.assertEqual(result.summary.n.tolist(), [2, 2])

    def test_nullable_numeric_missing_and_invalid_inputs(self):
        frame = pd.DataFrame({'score': pd.array([1, None, 3], dtype='Int64'), 'flag': [True, False, True], 'category': ['x'] * 3})
        masks = {'A': [True, True, False], 'B': [False, False, True]}
        self.assertEqual(numeric_columns(frame), ['score'])
        self.assertEqual(compare_metadata(frame, masks, ['score']).summary.value.tolist(), [1., 3.])
        for columns, chosen, statistic in [([], masks, 'mean'), (['category'], masks, 'mean'),
            (['score'], {'A': masks['A']}, 'mean'), (['score'], {'A': [False] * 3, 'B': masks['B']}, 'mean'),
            (['score'], {'A': [True], 'B': [True]}, 'mean'), (['score'], masks, 'sum')]:
            with self.assertRaises(ValueError):
                compare_metadata(frame, chosen, columns, statistic)
        with self.assertRaisesRegex(ValueError, 'No finite'):
            compare_metadata(pd.DataFrame({'v': [np.nan, np.inf]}), {'A': [True, False], 'B': [False, True]}, ['v'])

    def test_rest_excludes_all_defined_groups_and_produces_independent_pairs(self):
        frame = pd.DataFrame({'score': [0., 10., 20., 30., 40.]})
        defined = {'A': [True, False, False, False, False],
                   'B': [False, True, False, False, False],
                   'hidden': [False, False, True, False, False]}
        result = compare_metadata_to_rest(frame, {'A': defined['A'], 'B': defined['B']}, defined, ['score'])
        self.assertEqual(list(result), ['A', 'B'])
        for name, value in (('A', 0.), ('B', 10.)):
            self.assertEqual(result[name].groups, (name, 'The rest'))
            self.assertEqual(result[name].summary.value.tolist(), [value, 35.])
            self.assertEqual(result[name].summary.total.tolist(), [1, 2])
        # Clusters are comparison sets, not additional custom groups defining Rest.
        cluster = compare_metadata_to_rest(frame, {'Cluster 0': [False, False, False, True, False]}, defined, ['score'])
        self.assertEqual(cluster['Cluster 0'].summary.value.tolist(), [30., 35.])

    def test_empty_rest_no_custom_groups_and_name_collisions(self):
        frame = pd.DataFrame({'score': [1., 2., 3.]})
        with self.assertRaisesRegex(ValueError, 'The rest contains no observations'):
            compare_metadata_to_rest(frame, {'A': [True, False, False]}, {'All': [True] * 3}, ['score'])
        with self.assertRaisesRegex(ValueError, 'at least one'):
            compare_metadata_to_rest(frame, {}, {}, ['score'])
        with self.assertRaisesRegex(ValueError, 'align'):
            rest_mask({'A': [True]}, 3)
        self.assertTrue(rest_mask({}, 3).all())
        result = compare_metadata_to_rest(frame, {'The rest': [True, False, False]},
                                          {'The rest': [True, False, False]}, ['score'])
        self.assertEqual(result['The rest'].groups, ('The rest', 'The rest (outside defined groups)'))
        self.assertEqual(result['The rest'].summary.value.tolist(), [1., 2.5])

    def test_horizontal_plot_has_one_bar_series_per_group_and_fraction_axis(self):
        frame = pd.DataFrame({'T': [.1, .3, .8], 'B': [.9, .7, .2]})
        result = compare_metadata(frame, {'first': [True, True, False], 'second': [False, False, True]}, ['T', 'B'])
        chart = build_metadata_plot(result, fractions=True, unit='proportion')
        bars = [renderer for renderer in chart.renderers if renderer.glyph.__class__.__name__ == 'HBar']
        self.assertEqual(len(bars), 2)
        self.assertEqual(chart.y_range.factors, ['B', 'T'])
        self.assertEqual((chart.x_range.start, chart.x_range.end), (0, 1))
        self.assertEqual(chart.xaxis[0].axis_label, 'Mean proportion')
        self.assertTrue(any(isinstance(tool, HoverTool) for tool in chart.tools))
        self.assertEqual(bars[0].data_source.data['value'].tolist(), [.2, .8])
        document = Document()
        document.add_root(chart)
        document.to_json()


class MetadataTabTests(unittest.TestCase):
    def setup_comparison(self):
        layer = make_layer()
        names = layer.sample.adata.obs_names
        fractions = pd.DataFrame({'T cells': np.linspace(0, 1, len(names)),
                                  'B cells': np.linspace(1, 0, len(names))}, index=names)
        layer.deconvolution_manager.apply_result(DeconvolutionResult(fractions, ['g0'],
            metadata={'method': 'nnls', 'params': {}}), 'cell_type')
        groups = layer.groups_manager
        groups.add_mask_group('A', np.arange(len(names)) < 8, layer.roi_manager.doc)
        groups.add_mask_group('B', np.arange(len(names)) >= 8, layer.roi_manager.doc)
        manager = layer.metadata_comparison_manager
        manager.get_widgets(ImmediateDoc())
        manager._checkboxes['groups'].active = [0, 1]
        return layer, manager

    def run_comparison(self, manager):
        manager.run._trigger_event(ButtonClick(manager.run))
        wait_for(lambda: not manager.logger.running and not manager.run.disabled)
        self.assertIsNotNone(manager.result, manager.logger.div.text)

    def test_separate_tab_preserves_de_and_renders_results_below_map(self):
        layer, manager = self.setup_comparison()
        tabs = layer.get_widgets(ImmediateDoc())
        self.assertIn('DE', [tab.title for tab in tabs.tabs])
        self.assertIn('Metadata comparison', [tab.title for tab in tabs.tabs])
        manager._checkboxes['groups'].active = [0, 1]
        self.assertEqual(manager.preset.value, '1:fraction')
        self.assertEqual(manager.columns.value, ['deconv_1:T cells', 'deconv_1:B cells'])
        layer.manager.set_fetch_context(bbox=(0, 1, 0, 1), resolution_level='agg:2')
        self.run_comparison(manager)
        result = manager.result.summary.set_index(['group', 'column'])
        self.assertAlmostEqual(result.loc[('A', 'deconv_1:T cells'), 'value'], np.linspace(0, 1, 16)[:8].mean())
        self.assertEqual(result.loc[('A', 'deconv_1:T cells'), 'n'], 8)
        self.assertEqual(len(layer._figure_de_layout.children), 2)
        chart = layer._figure_de_layout.children[1]
        self.assertEqual(chart.y_range.factors, ['B cells', 'T cells'])
        self.assertIs(layer._figure_de_layout._comparison_owner, manager)
        # DE clearing does not erase a more recent metadata comparison.
        layer.de_manager._clear_de_results()
        self.assertEqual(len(layer._figure_de_layout.children), 2)
        layer.groups_manager.delete_group('A', manager.logger.doc)
        self.assertIsNone(manager.result)
        self.assertEqual(layer._figure_de_layout.children, [])

    def setup_rest_comparison(self):
        layer, manager = self.setup_comparison()
        groups = layer.groups_manager
        for name in list(groups._groups):
            groups.delete_group(name, manager.logger.doc)
        for name, start in [('A', 0), ('B', 4), ('Unselected', 8)]:
            groups.add_mask_group(name, (np.arange(16) >= start) & (np.arange(16) < start + 4), manager.logger.doc)
        groups._active_groups.discard('Unselected')
        groups._apply_visibility_mask()
        manager.include_rest.active = True
        return layer, manager

    def test_rest_uses_hidden_unselected_groups_and_one_plot_per_selected_set(self):
        layer, manager = self.setup_rest_comparison()
        manager._checkboxes['groups'].active = [0, 1]
        self.run_comparison(manager)
        self.assertEqual(list(manager.result), ['A', 'B'])
        self.assertEqual(len(layer._figure_de_layout.children), 3)
        titles = [layer._figure_de_layout.children[index].title.text for index in (1, 2)]
        self.assertEqual(titles, ['A vs The rest', 'B vs The rest'])
        for result in manager.result.values():
            rest = result.summary[result.summary.group == 'The rest']
            self.assertTrue((rest.total == 4).all())
            self.assertAlmostEqual(rest[rest.column == 'deconv_1:T cells'].value.iloc[0], np.linspace(0, 1, 16)[12:].mean())
        # New unselected groups change the remainder and invalidate saved comparisons.
        layer.groups_manager.add_mask_group('Extra', np.arange(16) == 12, manager.logger.doc)
        self.assertIsNone(manager.result)
        self.assertEqual(layer._figure_de_layout.children, [])

    def test_rest_allows_one_selected_cluster_or_group(self):
        layer, manager = self.setup_rest_comparison()
        manager._checkboxes['groups'].active = [0]
        self.run_comparison(manager)
        self.assertEqual(list(manager.result), ['A'])
        self.assertEqual(len(layer._figure_de_layout.children), 2)
        clustering = layer.clustering_manager
        clustering._cluster_columns['leiden'] = 'labels'
        clustering.cluster_col = 'labels'
        clustering._cluster_labels = ['0', '1']
        clustering._color_map = {'0': '#123456', '1': '#654321'}
        layer.sample.adata.obs['labels'] = ['0'] * 8 + ['1'] * 8
        manager.refresh_clusters()
        manager._checkboxes['groups'].active = []
        manager._checkboxes['clusters'].active = [1]
        self.run_comparison(manager)
        self.assertEqual(list(manager.result), ['Cluster 1'])
        self.assertEqual(layer._figure_de_layout.children[1].title.text, 'Cluster 1 vs The rest')
        document = Document()
        document.add_root(layer._figure_de_layout)
        document.to_json()

    def test_empty_rest_reports_error_and_releases_run_button(self):
        layer, manager = self.setup_comparison()  # A and B already cover every observation.
        manager.include_rest.active = True
        manager._checkboxes['groups'].active = [0]
        manager.run._trigger_event(ButtonClick(manager.run))
        wait_for(lambda: not manager.logger.running and not manager.run.disabled)
        self.assertIn('The rest contains no observations', manager.logger.div.text)
        self.assertIsNone(manager.result)
        self.assertEqual(layer._figure_de_layout.children, [])

    def test_new_deconvolution_runs_refresh_choices_without_recreating_layer(self):
        layer, manager = self.setup_comparison()
        fractions = layer.sample.adata.obsm['deconvolution'].drop(columns='most_abundant')
        layer.deconvolution_manager.apply_result(DeconvolutionResult(fractions, ['g0'], metadata={'method': 'music'}), 'cell_type')
        self.assertIn('2:fraction', dict(manager.preset.options))
        manager.preset.value = '2:fraction'
        self.assertEqual(manager.columns.value, ['deconv_2:T cells', 'deconv_2:B cells'])
        manager.preset.value = 'numeric'
        self.assertIn('deconv_1:T cells', dict(manager.columns.options))
        self.assertNotIn('deconv_1:most_abundant', dict(manager.columns.options))

    def test_three_sets_can_mix_groups_and_clusters_and_reclustering_invalidates(self):
        layer, manager = self.setup_comparison()
        clustering = layer.clustering_manager
        clustering._cluster_columns['leiden'] = 'labels'
        clustering.cluster_col = 'labels'
        clustering._cluster_labels = ['0', '1']
        clustering._color_map = {'0': '#123456', '1': '#654321'}
        layer.sample.adata.obs['labels'] = ['0'] * 8 + ['1'] * 8
        manager.refresh_clusters()
        manager._checkboxes['clusters'].active = [1]
        self.run_comparison(manager)
        self.assertEqual(manager.result.groups, ('A', 'B', 'Cluster 1'))
        means = manager.result.summary.set_index(['group', 'column'])
        self.assertEqual(means.loc[('B', 'deconv_1:T cells'), 'value'],
                         means.loc[('Cluster 1', 'deconv_1:T cells'), 'value'])
        self.assertEqual(len(layer._figure_de_layout.children[1].renderers), 3)
        # Invalidating metadata must not erase a subsequently displayed DE result.
        previous = Div(text='New DE result')
        layer._figure_de_layout.children = [previous]
        layer._figure_de_layout._comparison_owner = layer.de_manager
        manager.refresh_clusters()
        self.assertIsNone(manager.result)
        self.assertEqual(layer._figure_de_layout.children, [previous])

    def test_membership_changes_during_work_discard_results_and_release_button(self):
        layer, manager = self.setup_comparison()
        started, release = threading.Event(), threading.Event()
        def slow(*args):
            started.set()
            release.wait(5)
            return compare_metadata(*args)
        with patch('hd_viewer.board.layout.layer.metadata_comparison.compare_metadata', side_effect=slow):
            manager.run._trigger_event(ButtonClick(manager.run))
            self.assertTrue(started.wait(5))
            layer.sample.adata.obs['grp_A'] = np.arange(16) < 4
            layer.groups_manager._notify_group_change()
            release.set()
            wait_for(lambda: not manager.logger.running and not manager.run.disabled)
        self.assertIsNone(manager.result)
        self.assertIn('discarded', manager.logger.div.text)
        self.assertEqual(layer._figure_de_layout.children, [])

    def test_invalid_selection_and_failures_keep_existing_de_output(self):
        layer, manager = self.setup_comparison()
        previous = Div(text='Existing DE result')
        layer._figure_de_layout.children = [previous]
        layer._figure_de_layout._comparison_owner = layer.de_manager
        manager._checkboxes['groups'].active = [0]
        manager.run._trigger_event(ButtonClick(manager.run))
        self.assertIn('at least two', manager.logger.div.text)
        manager._checkboxes['groups'].active = [0, 1]
        with patch('hd_viewer.board.layout.layer.metadata_comparison.compare_metadata', side_effect=ValueError('Invalid metadata')):
            manager.run._trigger_event(ButtonClick(manager.run))
            wait_for(lambda: not manager.logger.running and not manager.run.disabled)
        self.assertIn('Invalid metadata', manager.logger.div.text)
        self.assertEqual(layer._figure_de_layout.children, [previous])
