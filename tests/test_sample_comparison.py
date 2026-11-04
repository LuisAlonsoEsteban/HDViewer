"""compare_samples honours the caller's sample order and per-sample keys; the panel labels results accordingly."""
import unittest
from types import SimpleNamespace
from unittest.mock import Mock

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from hd_viewer.board.operations.comparison import compare_samples, oriented_name
from hd_viewer.startup.sample._sample import Sample
from hd_viewer.help.comparison import LOCAL_HELP


def _grid(n, size, prefix):
    names, boxes = [], []
    for i in range(n):
        for j in range(n):
            names.append(f'{prefix}{i}_{j}')
            boxes.append(shapely.box(j * size, i * size, (j + 1) * size, (i + 1) * size))
    return names, gpd.GeoDataFrame({'name': names}, geometry=boxes, index=names)


def _samples():
    rng = np.random.default_rng(0)
    genes = [f'g{i}' for i in range(12)]
    # 3 x 3 spots of 10 units, each holding 2 x 2 bins of 5 units
    spot_names, spot_gdf = _grid(3, 10., 's')
    bin_names, bin_gdf = _grid(6, 5., 'b')
    bin_counts = rng.poisson(3, (len(bin_names), len(genes))).astype(float)
    spot_of_bin = [f's{int(i) // 2}_{int(j) // 2}' for i, j in (n[1:].split('_') for n in bin_names)]
    spot_counts = (pd.DataFrame(bin_counts, index=spot_of_bin).groupby(level=0).sum()
                   .reindex(spot_names).to_numpy() + rng.poisson(1, (len(spot_names), len(genes))))

    spots = ad.AnnData(spot_counts, obs=pd.DataFrame(index=spot_names), var=pd.DataFrame(index=genes))
    bins = ad.AnnData(bin_counts, obs=pd.DataFrame(index=bin_names), var=pd.DataFrame(index=genes))
    # Asymmetric keys: each one exists in one sample only
    spots.obs['spot_cluster'] = pd.Categorical(['a', 'b', 'a'] * 3)
    spots.obs['spot_type'] = ['T', 'B', 'T'] * 3
    bins.obs['bin_cluster'] = pd.Categorical(np.where(np.arange(len(bin_names)) % 3, 'a', 'b'))
    bins.obs['bin_type'] = np.where(np.arange(len(bin_names)) % 2, 'T', 'B')
    t = np.linspace(0, 1, len(spot_names))
    spots.obsm['spot_deconv'] = pd.DataFrame({'T': t, 'B': 1 - t, 'most_abundant': 'T'}, index=spot_names)
    u = np.linspace(0, 1, len(bin_names))
    bins.obsm['bin_deconv'] = pd.DataFrame({'T': u, 'B': 1 - u}, index=bin_names)
    return Sample(spot_gdf, spots, on='name'), Sample(bin_gdf, bins, on='name')


class CompareSamplesOrderTests(unittest.TestCase):
    def setUp(self):
        self.spots, self.bins = _samples()
        self.kwargs = dict(min_coverage=.5, n_hvg=12, k_neighbors=2)

    def test_reversed_order_swaps_keys_with_samples_and_labels_results_in_caller_order(self):
        canonical = compare_samples(self.spots, self.bins, **self.kwargs,
                                    sample1_cluster_key='spot_cluster', sample2_cluster_key='bin_cluster',
                                    sample1_celltype_key='spot_deconv', sample2_celltype_key='bin_type')
        # Sample 1 = bins, Sample 2 = spots: previously KeyError ('bin_cluster' looked up in the spots)
        reversed_ = compare_samples(self.bins, self.spots, **self.kwargs,
                                    sample1_cluster_key='bin_cluster', sample2_cluster_key='spot_cluster',
                                    sample1_celltype_key='bin_type', sample2_celltype_key='spot_deconv')
        self.assertEqual(canonical['reference_sample'], 1)
        self.assertEqual(reversed_['reference_sample'], 2)

        a, b = canonical['local'], reversed_['local']
        self.assertEqual(len(a), 9)
        pd.testing.assert_index_equal(a.index, b.index)
        # Per-sample columns follow the caller's order: sample1 = bins in the reversed call
        np.testing.assert_allclose(b['total_counts_sample1'], a['total_counts_sample2'])
        np.testing.assert_allclose(b['total_counts_sample2'], a['total_counts_sample1'])
        self.assertEqual(b['cluster_sample2'].tolist(), a['cluster_sample1'].tolist())
        self.assertEqual(b['cluster_sample1_majority'].tolist(), a['cluster_sample2_majority'].tolist())
        self.assertEqual(b['n_sample1'].tolist(), a['n_sample2'].tolist())
        np.testing.assert_allclose(b['celltype_jensen_shannon'], a['celltype_jensen_shannon'])
        self.assertTrue(np.isfinite(a['celltype_jensen_shannon']).all())
        self.assertNotIn('cluster_sample1', b.columns)

        ga, gb = canonical['global'], reversed_['global']
        self.assertEqual(gb['n_sample1_matched'], ga['n_sample2_matched'])
        self.assertEqual(gb['n_sample2'], ga['n_sample1'])
        self.assertEqual(gb['tissue_area_sample1'], ga['tissue_area_sample2'])
        self.assertEqual(gb['cluster_ari'], ga['cluster_ari'])
        self.assertEqual(gb['cluster_sample1_purity_mean'], ga['cluster_sample2_purity_mean'])
        np.testing.assert_allclose(reversed_['genes']['morans_i_sample2'], canonical['genes']['morans_i_sample1'])
        np.testing.assert_allclose(reversed_['genes']['zero_fraction_sample1'], canonical['genes']['zero_fraction_sample2'])
        # Canonical names (used by the help texts) are recovered from either orientation
        self.assertEqual(oriented_name('n_sample1_matched', True), 'n_sample2_matched')
        self.assertEqual(oriented_name('n_sample1_matched', False), 'n_sample1_matched')
        self.assertEqual({oriented_name(c, True) for c in b.columns}, set(a.columns))
        self.assertEqual({oriented_name(k, True) for k in gb}, set(ga))

    def test_deconvolution_on_high_res_and_labels_on_low_res(self):
        """Sample 1 = bins with a deconvolution obsm key, Sample 2 = spots with labels."""
        result = compare_samples(self.bins, self.spots, **self.kwargs,
                                 sample1_celltype_key='bin_deconv', sample2_celltype_key='spot_type')
        js = result['local']['celltype_jensen_shannon']
        self.assertEqual(len(js), 9)
        self.assertTrue(np.isfinite(js).all())
        self.assertEqual(set(result['global']['celltype_spatial_pearson']), {'T', 'B'})

    def test_unknown_celltype_key_reports_which_sample(self):
        with self.assertRaisesRegex(ValueError, 'high-res'):
            compare_samples(self.bins, self.spots, **self.kwargs,
                            sample1_celltype_key='spot_type', sample2_celltype_key='spot_deconv')


class ComparisonPanelOrderTests(unittest.TestCase):
    def _component(self, name, alias, sample):
        class Component:  # hashable stand-in for a layer component
            pass
        component = Component()
        component.__dict__.update(name=name, alias=alias, sample=sample, autocomplete_input=None,
                                  manager=SimpleNamespace(cache=pd.DataFrame()),
                                  refresh_completions=Mock(), get_colnames=lambda: [])
        return component

    def test_panel_options_follow_user_order_and_results_land_on_low_res_sample(self):
        from bokeh.document import Document
        from hd_viewer.board.layout.comparison import ComparisonPanel

        spots, bins = _samples()
        panel = ComparisonPanel(Document())
        spot_c, bin_c = self._component('spots', 'Spot', spots), self._component('bins', 'Bin', bins)
        panel.set_components([spot_c, bin_c])
        panel.sample1_select.value, panel.sample2_select.value = 'bins', 'spots'
        values = lambda select: [o[0] if isinstance(o, tuple) else o for o in select.options]
        self.assertIn('bin_cluster', panel.sample1_cluster_select.options)
        self.assertIn('bin_deconv', values(panel.sample1_celltype_select))
        self.assertIn('spot_type', values(panel.sample2_celltype_select))
        self.assertIn('spot_deconv', values(panel.sample2_celltype_select))

        result = compare_samples(bins, spots, min_coverage=.5, n_hvg=12, k_neighbors=2,
                                 sample1_cluster_key='bin_cluster', sample2_cluster_key='spot_cluster')
        panel._on_result(bin_c, spot_c, result)
        self.assertIs(panel.sample1_component, spot_c)  # local metrics belong to the spots
        self.assertTrue(panel.swapped)
        obs = spots.adata.obs
        self.assertIn('cmp_total_counts_sample2', obs)
        self.assertNotIn('cmp_cluster_sample2', obs)  # the spots' own clusters are not duplicated
        self.assertIn('cmp_cluster_sample1_majority', obs)
        self.assertNotIn('cmp_total_counts_sample1', bins.adata.obs)
        # Gene table keeps low/high titles: its "low" columns hold the spots' values
        np.testing.assert_allclose(panel.gene_source.data['morans_i_sample1'],
                                   result['genes'].sort_values('spatial_pearson', ascending=False)['morans_i_sample2'])
        # Help is looked up with the canonical (low-res = sample1) name
        data = panel.global_source.data
        i = data['metric'].index('n_sample2')
        self.assertIn('low-res', data['help'][i])
        panel.metric_select.value = 'total_counts_sample2'
        self.assertEqual(panel.metric_help.content, LOCAL_HELP['total_counts_sample1'])
        self.assertIn('Compared 9 spots with 36 bins', panel.status.div.text)


if __name__ == '__main__':
    unittest.main()
