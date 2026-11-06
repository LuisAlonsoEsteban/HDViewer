"""Descriptive comparison of numeric metadata between groups/clusters or against The rest."""
import json
import threading
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.models import HoverTool

from hd_viewer.analysis.deconvolution import DeconvolutionResult
from hd_viewer.analysis.metadata_comparison import (build_metadata_plot, compare_metadata,
                                                    compare_metadata_to_rest, numeric_columns, rest_mask)
from hd_viewer.api.routers import metadata as metadata_router
from hd_viewer.core.deconvolution import publish_result
from tests.test_roi_manager import Api

T_CELLS = np.linspace(0, 1, 16)
COLUMNS = ['deconv_1:T cells', 'deconv_1:B cells']


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
        # Clusters are comparison sets, not additional custom groups defining Rest; like DE vs Rest,
        # a cluster's own observations are never part of its rest.
        cluster = compare_metadata_to_rest(frame, {'Cluster 0': [False, False, False, True, False]}, defined, ['score'])
        self.assertEqual(cluster['Cluster 0'].summary.value.tolist(), [30., 40.])
        self.assertEqual(cluster['Cluster 0'].summary.total.tolist(), [1, 1])

    def test_rest_excludes_the_compared_cluster_even_when_it_overlaps_groups(self):
        frame = pd.DataFrame({'score': [0., 10., 20., 30., 40., 50.]})
        defined = {'A': [True, True, False, False, False, False]}
        # Cluster 3 straddles group A and the ungrouped observations
        cluster = [False, True, True, True, False, False]
        result = compare_metadata_to_rest(frame, {'Cluster 3': cluster, 'A': defined['A']}, defined, ['score'])
        rest = result['Cluster 3'].summary.set_index('group').loc['The rest']
        self.assertEqual((rest.value, rest.total), (45., 2))
        self.assertEqual(result['Cluster 3'].summary.value.iloc[0], 20.)
        # Each set has its own rest: A's rest still contains the cluster's ungrouped observations
        rest = result['A'].summary.set_index('group').loc['The rest']
        self.assertEqual((rest.value, rest.total), (35., 4))

    def test_rest_empty_once_the_compared_cluster_is_removed(self):
        frame = pd.DataFrame({'score': [1., 2., 3.]})
        defined = {'A': [True, False, False]}
        with self.assertRaisesRegex(ValueError, "The rest contains no observations outside the defined groups and 'Cluster 0'"):
            compare_metadata_to_rest(frame, {'Cluster 0': [False, True, True]}, defined, ['score'])
        with self.assertRaisesRegex(ValueError, 'align'):
            compare_metadata_to_rest(frame, {'Cluster 0': [False, True]}, defined, ['score'])

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


class MetadataApiTests(unittest.TestCase):
    """The 16-bin grid with one deconvolution run (T cells rising from 0 to 1 over the bins)."""
    def setup_comparison(self):
        api = Api(self)
        sample = api.view.sample
        names = sample.adata.obs_names
        fractions = pd.DataFrame({'T cells': T_CELLS, 'B cells': T_CELLS[::-1]}, index=names)
        publish_result(sample, DeconvolutionResult(fractions, ['g0'], metadata={'method': 'nnls', 'params': {}}),
                       'cell_type')
        self.add_group(api, 'A', np.arange(16) < 8)
        self.add_group(api, 'B', np.arange(16) >= 8)
        return api

    @staticmethod
    def add_group(api, name, mask):
        api.view.groups.add_fixed_group(name, mask)
        api.view.drop_comparisons(groups_changed=True)     # as the API does on every new group

    @staticmethod
    def with_clusters(api, labels):
        clustering = api.view.clustering
        api.view.sample.adata.obs['labels'] = labels
        clustering._cluster_columns['leiden'] = clustering.cluster_col = 'labels'
        clustering._build_cluster_state()
        clustering._color_map.update({'0': '#123456', '1': '#654321'})

    def compare(self, api, **body):
        """POST a comparison; returns (response, comparisons computed by the router)."""
        computed = []

        def recorded(function):
            def wrapper(*args, **kwargs):
                result = function(*args, **kwargs)
                computed.append(result)
                return result
            return wrapper
        body.setdefault('preset', '1:fraction')
        body.setdefault('columns', COLUMNS)
        with patch.object(metadata_router, 'compare_metadata', recorded(compare_metadata)), \
                patch.object(metadata_router, 'compare_metadata_to_rest', recorded(compare_metadata_to_rest)):
            response = api.post('metadata/compare', json=body)
        return response, (computed[0] if computed else None)

    def ready(self, api):
        return api.ok(api.get('metadata/options'))['ready']

    def plots(self, api):
        return json.dumps(api.ok(api.get('metadata/plots')))

    def test_summaries_use_all_native_observations_and_follow_group_changes(self):
        api = self.setup_comparison()
        options = api.ok(api.get('metadata/options'))
        self.assertEqual(options['default_preset'], '1:fraction')
        preset = next(p for p in options['presets'] if p['key'] == '1:fraction')
        self.assertEqual([c['name'] for c in preset['columns']], COLUMNS)
        self.assertEqual([c['label'] for c in preset['columns']], ['T cells', 'B cells'])
        # Independent of the displayed aggregation level and viewport
        api.view.rows.set_fetch_context(bbox=(0, 1, 0, 1), resolution_level='agg:2')
        api.view.rows.fetch_data()
        response, result = self.compare(api, groups=['A', 'B'])
        self.assertEqual(api.ok(response), dict(ready=True, plots=1))
        summary = result.summary.set_index(['group', 'column'])
        self.assertAlmostEqual(summary.loc[('A', 'deconv_1:T cells'), 'value'], T_CELLS[:8].mean())
        self.assertEqual(summary.loc[('A', 'deconv_1:T cells'), 'n'], 8)
        self.assertIn('"B cells"', self.plots(api))                       # cell-type labels on the axis
        # A DE result being dropped (re-clustering) does not erase this comparison of groups
        api.view.drop_comparisons(clusters=True)
        self.assertTrue(self.ready(api))
        api.ok(api.delete('groups/A'))
        self.assertFalse(self.ready(api))
        self.assertEqual(api.get('metadata/plots').status_code, 404)

    def setup_rest_comparison(self):
        api = self.setup_comparison()
        for name in list(api.view.groups._groups):
            api.ok(api.delete(f'groups/{name}'))
        for name, start in [('A', 0), ('B', 4), ('Unselected', 8)]:
            self.add_group(api, name, (np.arange(16) >= start) & (np.arange(16) < start + 4))
        self.assertEqual(api.patch('groups/Unselected', json=dict(visible=False)).status_code, 200)
        return api

    def test_rest_uses_hidden_unselected_groups_and_one_plot_per_selected_set(self):
        api = self.setup_rest_comparison()
        response, result = self.compare(api, groups=['A', 'B'], include_rest=True)
        self.assertEqual(api.ok(response)['plots'], 2)
        self.assertEqual(list(result), ['A', 'B'])
        plots = self.plots(api)
        self.assertIn('A vs The rest', plots)
        self.assertIn('B vs The rest', plots)
        for comparison in result.values():
            rest = comparison.summary[comparison.summary.group == 'The rest']
            self.assertTrue((rest.total == 4).all())
            self.assertAlmostEqual(rest[rest.column == 'deconv_1:T cells'].value.iloc[0], T_CELLS[12:].mean())
        # New unselected groups change the remainder and invalidate saved comparisons.
        roi = api.roi(kind='box', box=[.2, 3.2, .8, 3.8])                  # bin b0_3 (12)
        api.roi_group(roi)
        self.assertFalse(self.ready(api))

    def test_rest_allows_one_selected_cluster_or_group(self):
        api = self.setup_rest_comparison()
        response, result = self.compare(api, groups=['A'], include_rest=True)
        self.assertEqual(api.ok(response)['plots'], 1)
        self.assertEqual(list(result), ['A'])
        # Cluster 1 must not cover every ungrouped observation (12-15), or its rest would be empty
        self.with_clusters(api, ['1'] * 8 + ['0'] * 8)
        response, result = self.compare(api, clusters=['Cluster 1'], include_rest=True)
        self.assertEqual(api.ok(response)['plots'], 1)
        self.assertEqual(list(result), ['Cluster 1'])
        self.assertIn('Cluster 1 vs The rest', self.plots(api))
        self.assertEqual(api.view.metadata['compared']['chosen'], [('clusters', 'Cluster 1')])

    def test_cluster_vs_rest_excludes_the_cluster(self):
        api = self.setup_rest_comparison()                 # groups cover 0-11; 12-15 are ungrouped
        self.with_clusters(api, ['0'] * 14 + ['1'] * 2)
        response, result = self.compare(api, clusters=['Cluster 1'], include_rest=True)
        api.ok(response)
        summary = result['Cluster 1'].summary
        rest = summary[(summary.group == 'The rest') & (summary.column == 'deconv_1:T cells')].iloc[0]
        self.assertEqual(rest.total, 2)
        self.assertAlmostEqual(rest.value, T_CELLS[12:14].mean())
        # A cluster covering all ungrouped observations leaves an empty rest: a clear error
        api.view.drop_comparisons(clusters=True)                  # as re-clustering does
        self.assertFalse(self.ready(api))
        self.with_clusters(api, ['0'] * 12 + ['1'] * 4)
        response, _ = self.compare(api, clusters=['Cluster 1'], include_rest=True)
        self.assertEqual(response.status_code, 422)
        self.assertIn('The rest contains no observations', response.json()['detail'])
        self.assertIn('Cluster 1', response.json()['detail'])
        self.assertFalse(self.ready(api))

    def test_empty_rest_reports_error(self):
        api = self.setup_comparison()                      # A and B already cover every observation.
        response, _ = self.compare(api, groups=['A'], include_rest=True)
        self.assertEqual(response.status_code, 422)
        self.assertIn('The rest contains no observations', response.json()['detail'])
        self.assertFalse(self.ready(api))

    def test_new_deconvolution_runs_refresh_choices(self):
        api = self.setup_comparison()
        sample = api.view.sample
        fractions = sample.adata.obsm['deconvolution'].drop(columns='most_abundant')
        publish_result(sample, DeconvolutionResult(fractions, ['g0'], metadata={'method': 'music'}), 'cell_type')
        options = api.ok(api.get('metadata/options'))
        presets = {p['key']: p for p in options['presets']}
        self.assertIn('2:fraction', presets)
        self.assertEqual(options['default_preset'], '2:fraction')
        self.assertEqual([c['name'] for c in presets['2:fraction']['columns']],
                         ['deconv_2:T cells', 'deconv_2:B cells'])
        numeric = [c['name'] for c in presets['numeric']['columns']]
        self.assertIn('deconv_1:T cells', numeric)
        self.assertNotIn('deconv_1:most_abundant', numeric)

    def test_three_sets_can_mix_groups_and_clusters_and_reclustering_invalidates(self):
        api = self.setup_comparison()
        self.with_clusters(api, ['0'] * 8 + ['1'] * 8)
        response, result = self.compare(api, groups=['A', 'B'], clusters=['Cluster 1'])
        api.ok(response)
        self.assertEqual(result.groups, ('A', 'B', 'Cluster 1'))
        means = result.summary.set_index(['group', 'column'])
        self.assertEqual(means.loc[('B', 'deconv_1:T cells'), 'value'],
                         means.loc[('Cluster 1', 'deconv_1:T cells'), 'value'])
        self.assertIn('#654321', self.plots(api))                      # Cluster 1's colour
        # Re-clustering invalidates the comparison, but not a DE result of groups only.
        api.view.de.update(header='DE vs Rest', merged=pd.DataFrame(), compared=('groups', ['A']))
        api.view.drop_comparisons(clusters=True)
        self.assertFalse(self.ready(api))
        self.assertEqual(api.view.de['header'], 'DE vs Rest')

    def test_membership_changes_during_work_discard_results(self):
        api = self.setup_comparison()
        started, release = threading.Event(), threading.Event()

        def slow(*args):
            started.set()
            release.wait(5)
            return compare_metadata(*args)
        responses = []
        with patch.object(metadata_router, 'compare_metadata', side_effect=slow):
            worker = threading.Thread(target=lambda: responses.append(api.post('metadata/compare', json=dict(
                groups=['A', 'B'], preset='1:fraction', columns=COLUMNS))))
            worker.start()
            self.assertTrue(started.wait(5))
            api.ok(api.delete('groups/A'))
            release.set()
            worker.join(30)
        self.assertFalse(self.ready(api))
        self.assertEqual(responses[0].status_code, 409)

    def test_invalid_selection_and_failures_keep_existing_results(self):
        api = self.setup_comparison()
        api.view.de.update(header='Existing DE result', compared=('groups', ['A']))
        response, _ = self.compare(api, groups=['A'])
        self.assertEqual(response.status_code, 422)
        self.assertIn('at least two', response.json()['detail'])
        api.ok(self.compare(api, groups=['A', 'B'])[0])
        with patch.object(metadata_router, 'compare_metadata', side_effect=ValueError('Invalid metadata')):
            response = api.post('metadata/compare', json=dict(groups=['A', 'B'], preset='1:fraction', columns=COLUMNS))
        self.assertEqual(response.status_code, 422)
        self.assertIn('Invalid metadata', response.json()['detail'])
        self.assertEqual(api.view.de['header'], 'Existing DE result')
        self.assertTrue(self.ready(api))                               # the previous comparison stays
        for body in (dict(statistic='sum'), dict(preset='nope'), dict(columns=['g0'])):
            with self.subTest(body=body):
                self.assertEqual(self.compare(api, groups=['A', 'B'], **body)[0].status_code, 422)


if __name__ == '__main__':
    unittest.main()
