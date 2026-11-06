"""compare_samples honours the caller's sample order and per-sample keys; results are labelled accordingly."""
import unittest
import time
from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from fastapi.testclient import TestClient

from hd_viewer.analysis.comparison import compare_samples, oriented_name
from hd_viewer.startup.sample._sample import Sample
from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.api.views import build_views
from hd_viewer.core import comparison as core
from hd_viewer.help.comparison import LOCAL_HELP
from hd_viewer.settings.config import comparison_kwargs


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


class ComparisonOrderTests(unittest.TestCase):
    def test_options_follow_each_sample_and_results_land_on_the_low_res_sample(self):
        spots, bins = _samples()
        self.assertIn('bin_cluster', core.label_columns(bins))
        self.assertNotIn('spot_cluster', core.label_columns(bins))
        values = lambda sample: [value for value, _ in core.celltype_options(sample)]
        self.assertIn('bin_deconv', values(bins))
        self.assertIn('spot_type', values(spots))
        self.assertEqual(core.celltype_options(spots)[0], ('spot_deconv', 'spot_deconv (proportions)'))
        self.assertEqual(core.by_area({'bins': bins, 'spots': spots}), ['spots', 'bins'])

        # Sample 1 = bins, Sample 2 = spots: the spots are still the low-res reference
        result = compare_samples(bins, spots, min_coverage=.5, n_hvg=12, k_neighbors=2,
                                 sample1_cluster_key='bin_cluster', sample2_cluster_key='spot_cluster')
        swapped = result['reference_sample'] == 2
        self.assertTrue(swapped)
        metrics = core.write_obs_columns(spots, result['local'], swapped)
        obs = spots.adata.obs
        self.assertIn('cmp_total_counts_sample2', obs)
        self.assertIn('total_counts_sample2', metrics)
        self.assertNotIn('cmp_cluster_sample2', obs)  # the spots' own clusters are not duplicated
        self.assertIn('cmp_cluster_sample1_majority', obs)
        self.assertNotIn('cmp_total_counts_sample1', bins.adata.obs)
        # Comparison metrics are not offered as cluster labels
        self.assertEqual(core.label_columns(spots), ['spot_cluster', 'spot_type'])
        # Help is looked up with the canonical (low-res = sample1) name
        rows = {row['metric']: row for row in core.flatten_global(result['global'], swapped)}
        self.assertIn('low-res', rows['n_sample2']['help'])
        self.assertEqual(core.local_help('total_counts_sample2', swapped), LOCAL_HELP['total_counts_sample1'])
        core.clear_obs_columns(spots)
        self.assertFalse([c for c in obs.columns if c.startswith(core.PREFIX)])


class ComparisonApiOrderTests(unittest.TestCase):
    def test_gene_table_and_help_are_oriented_low_res_first(self):
        spots, bins = _samples()
        registry = SessionRegistry(ttl=3600)
        client = TestClient(create_app(registry, sweep_seconds=3600))
        client.__enter__()
        self.addCleanup(client.__exit__, None, None, None)
        client.get('/api/load')
        session = registry.get(client.cookies[COOKIE])
        session.loaded = SimpleNamespace(samples={'bins': ('spots', bins, 'Bins'), 'spots': ('spots', spots, 'Spots')})
        session.views = build_views(session.loaded)
        started = client.post('/api/comparison', json=dict(sample1='bins', sample2='spots', min_coverage=.5, predicate='contains',
                                                           sample1_cluster_key='bin_cluster',
                                                           sample2_cluster_key='spot_cluster'))
        self.assertEqual(started.status_code, 202, started.text)
        for _ in range(1200):
            job = client.get('/api/comparison').json()['job']
            if job['state'] != 'running':
                break
            time.sleep(.05)
        self.assertEqual(job['state'], 'done', job)
        self.assertIn('Compared 9 spots with 36 spots', job['message'])
        result = client.get('/api/comparison/result').json()
        self.assertEqual((result['low'], result['high']), ('spots', 'bins'))
        # Gene table keeps low/high titles: its "low" columns hold the spots' values
        expected = compare_samples(bins, spots, **{k: v for k, v in comparison_kwargs().items()
                                                   if k not in ('min_coverage', 'predicate')},
                                   min_coverage=.5, predicate='contains', sample1_cluster_key='bin_cluster',
                                   sample2_cluster_key='spot_cluster')['genes']
        expected = expected.sort_values('spatial_pearson', ascending=False)
        keys = [column['key'] for column in result['gene_columns']]
        low = keys.index('morans_i_sample1') + 1
        self.assertEqual([row[0] for row in result['genes']], [str(g) for g in expected.index])
        np.testing.assert_allclose([row[low] for row in result['genes']], expected['morans_i_sample2'])
        self.assertEqual(result['metric_help']['total_counts_sample2'], LOCAL_HELP['total_counts_sample1'])
        metrics = {row['metric']: row for row in result['global_metrics']}
        self.assertIn('low-res', metrics['n_sample2']['help'])


if __name__ == '__main__':
    unittest.main()
