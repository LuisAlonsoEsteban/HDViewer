"""Algorithm dispatch, cache invalidation, native projection and the API contract."""
import time
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import adjusted_rand_score

from hd_viewer.core.clustering import SUPPORTED_ALGORITHMS, ClusteringState
from hd_viewer.startup.data.image.transform import SpatialTransform
from hd_viewer.startup.operations.bins2spots import Bins2Spots
from hd_viewer.startup.sample import Bins
from hd_viewer.startup.sample._sample import Sample
from tests.test_clustering_markers import api_client, make_layer


def prepared_sample():
    rng = np.random.default_rng(9)
    values = np.r_[rng.normal(-5, .12, (20, 3)), rng.normal(5, .12, (20, 3))]
    names = [f'c{i}' for i in range(len(values))]
    sample = Sample(pd.DataFrame({'index': names}),
                    ad.AnnData(np.exp(values), obs=pd.DataFrame(index=names)))
    sample.adata_processed = sample.adata.copy()
    sample.adata_processed.obsm['X_pca'] = values
    sample._step_params['pca'] = 3
    return sample


def aligned_bins(n=6):
    """A 6 x 6 grid of 2-um bins registered on an image (from test_spatial_alignment)."""
    theta = np.deg2rad(23.)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    basis = np.array([[4., 0.], [.4, 3.5]]) @ rotation.T
    rows, cols = np.indices((n, n))
    xy = np.column_stack((cols.ravel(), rows.ravel())) @ basis + [65, 70]
    names = [f'b{i}' for i in range(len(xy))]
    positions = pd.DataFrame(dict(array_col=cols.ravel(), array_row=rows.ravel(),
                                  pxl_col_in_fullres=xy[:, 0], pxl_row_in_fullres=xy[:, 1], in_tissue=1,
                                  index=names), index=names)
    transform = SpatialTransform.fit(positions, (140, 180, 3), 1.)
    counts = np.column_stack((np.arange(len(xy)) + 1, np.ones(len(xy))))
    data = ad.AnnData(counts, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['gene', 'one']))
    b2s = Bins2Spots(data, positions, pixels_per_bin=4, bins_size=2, spatial_transform=transform)
    return Bins(b2s.bin_squares, data, pixels_per_bin=4, bin_size=2, spatial_transform=transform)


class Counted:
    """Counts the fits of an estimator class: patch.object(cls, 'fit_predict', counted.method)."""
    def __init__(self, estimator):
        self.calls, original = 0, estimator.fit_predict

        def method(estimator, *args, **kwargs):
            self.calls += 1
            return original(estimator, *args, **kwargs)
        self.method = method


class ClusteringAlgorithmTests(unittest.TestCase):
    def test_all_algorithms_find_separated_populations_without_mutating_counts(self):
        for algorithm in SUPPORTED_ALGORITHMS:
            with self.subTest(algorithm=algorithm):
                sample = prepared_sample()
                state = ClusteringState(sample)
                before = sample.adata.X.copy()
                column = state.run_clustering(algorithm=algorithm, n_pcs=3, n_neighbors=8,
                    resolution=.1, n_clusters=2, eps=1, min_samples=3, min_cluster_size=5)
                labels = sample.adata.obs[column]
                self.assertGreater(adjusted_rand_score([0]*20 + [1]*20, labels), .95)
                np.testing.assert_array_equal(sample.adata.X, before)
                self.assertEqual(sample.adata.uns['clustering']['algorithm'], algorithm)
                self.assertEqual(state.cluster_col, column)

    def test_density_noise_and_marker_edge_cases(self):
        sample = prepared_sample()
        state = ClusteringState(sample)
        sample.adata_processed.obsm['X_pca'][-1] = [100, 100, 100]
        column = state.run_clustering(algorithm='dbscan', n_pcs=3, eps=1, min_samples=3)
        self.assertEqual(sample.adata.obs[column].iloc[-1], 'Noise')
        state.compute_markers()
        self.assertNotIn('Noise', state.get_marker_genes())
        self.assertEqual(len(state.get_marker_genes()), 2)
        state._build_cluster_state()
        self.assertEqual(state._result_summary(), '2 clusters found, 1 noise observations')
        for algorithm, options in [('dbscan', {'eps': 1e-10}), ('kmeans', {'n_clusters': 1})]:
            with self.subTest(algorithm=algorithm):
                state.run_clustering(algorithm=algorithm, n_pcs=3, **options)
                state.compute_markers()
                self.assertEqual(state.get_marker_genes(), {})

    def test_cache_tracks_method_parameters_and_marker_freshness(self):
        sample = prepared_sample()
        state = ClusteringState(sample)
        kmeans, agglomerative = Counted(KMeans), Counted(AgglomerativeClustering)
        with patch.object(KMeans, 'fit_predict', kmeans.method), \
                patch.object(AgglomerativeClustering, 'fit_predict', agglomerative.method):
            state.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=2)
            state.compute_markers()
            self.assertTrue(state.get_marker_genes())
            state.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=2, eps=4)   # eps: not a kmeans option
            self.assertEqual(kmeans.calls, 1)
            self.assertTrue(state.get_marker_genes())
            state.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=3)
            self.assertEqual(kmeans.calls, 2)
            self.assertEqual(state.get_marker_genes(), {})
            state.run_clustering(algorithm='agglomerative', n_pcs=3, n_clusters=3)
            self.assertEqual(agglomerative.calls, 1)
        self.assertEqual(kmeans.calls, 2)

    def test_aggregate_algorithms_project_identical_labels_to_native_members(self):
        bins = aligned_bins()
        state = ClusteringState(bins)
        aggregate = bins._get_aggregated_sample(2)
        aggregate.adata_processed = aggregate.adata.copy()
        aggregate.adata_processed.obsm['X_pca'] = np.arange(aggregate.adata.n_obs*2).reshape(-1, 2)
        aggregate._step_params['pca'] = 2
        column = state.run_clustering(algorithm='kmeans', n_clusters=2, n_pcs=2,
                    resolution_level=bins.get_resolution_level_for_factor(2), compute_markers=False)
        projected = bins.aggregation_names(2).map(aggregate.adata.obs['kmeans_2'].astype(str))
        self.assertEqual(bins.adata.obs[column].tolist(), projected.tolist())
        self.assertEqual(column, 'kmeans_2')
        self.assertEqual(bins._clustering_factor, 2)
        self.assertEqual(state.get_cluster_factor(column), 2)
        self.assertIs(state._get_active_analysis_sample(), aggregate)

    def test_noise_visibility_and_hover(self):
        view = make_layer()
        clustering = view.clustering
        obs = view.sample.adata.obs
        obs['dbscan_1'] = ['Noise'] + ['0'] * (view.sample.adata.n_obs - 1)
        clustering._cluster_columns['dbscan'] = clustering.cluster_col = 'dbscan_1'
        clustering._build_cluster_state()
        self.assertEqual(clustering._cluster_labels, ['0', 'Noise'])
        clustering._active_clusters = {'Noise'}
        self.assertEqual(clustering.visible_mask(obs['dbscan_1']).sum(), 1)
        # A row without a cluster label (left out of the clustering) stays visible
        self.assertEqual(clustering.visible_mask([None, 'Noise', '0']).tolist(), [True, True, False])

        view.update_viewport(-1, 181, -1, 2, 900, 900)
        view.colour_by_categories('dbscan_1', clustering._color_map, 'clusters')
        self.assertEqual(int(view.visibility().sum()), 1)
        self.assertIn(['dbscan_1', '@{dbscan_1}'], view.tooltips())


class ClusteringApiTests(unittest.TestCase):
    def test_selected_algorithm_and_its_options_reach_the_run(self):
        client, session = api_client(self, make_layer().sample)
        clustering = session.views['bins'].clustering
        with patch.object(clustering, 'compute', side_effect=RuntimeError('stop')) as compute:
            response = client.post('/api/layers/bins/clustering',
                                   json=dict(algorithm='hdbscan', min_cluster_size=12, resolution=.5))
            self.assertEqual(response.status_code, 202)
            for _ in range(600):
                if client.get('/api/layers/bins/clusters').json()['job']['state'] != 'running':
                    break
                time.sleep(.05)
        self.assertEqual(compute.call_args.args[0], .5)
        self.assertEqual(compute.call_args.kwargs['algorithm'], 'hdbscan')
        self.assertEqual(compute.call_args.kwargs['min_cluster_size'], 12)
        self.assertNotIn('eps', compute.call_args.kwargs)        # unset options keep their defaults
        job = client.get('/api/layers/bins/clusters').json()['job']
        self.assertEqual((job['state'], job['message']), ('failed', 'stop'))
        self.assertEqual(len(client.get('/api/clustering/options').json()['algorithms']), 6)

    def test_invalid_algorithm_and_cluster_counts(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported'):
            ClusteringState(prepared_sample()).run_clustering(algorithm='invalid')
        client, _ = api_client(self, make_layer().sample)
        self.assertEqual(client.post('/api/layers/bins/clustering', json=dict(algorithm='invalid')).status_code, 422)
        for count in (0, 50):
            with self.subTest(count=count), self.assertRaises(ValueError):
                ClusteringState(prepared_sample()).run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=count)


if __name__ == '__main__':
    unittest.main()
