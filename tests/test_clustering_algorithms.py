"""Algorithm dispatch, cache invalidation, native projection and UI contracts."""
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.models import Select
from sklearn.metrics import adjusted_rand_score

from hd_viewer.startup.sample._sample import Sample
try:
    from src.utils.clustering import fit_clusters  # removed in refactor (f4265e9)
except ImportError:
    fit_clusters = None
from tests.test_clustering_markers import make_layer
from tests.test_spatial_alignment import fixture


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


@unittest.skipIf(fit_clusters is None, "fit_clusters / Sample.run_clustering removed in refactor")
class ClusteringAlgorithmTests(unittest.TestCase):
    def test_all_algorithms_find_separated_populations_without_mutating_counts(self):
        for algorithm in ('leiden', 'kmeans', 'dbscan', 'hdbscan', 'gaussian_mixture', 'agglomerative'):
            with self.subTest(algorithm=algorithm):
                sample = prepared_sample()
                before = sample.adata.X.copy()
                column = sample.run_clustering(algorithm=algorithm, n_pcs=3, n_neighbors=8,
                    resolution=.1, n_clusters=2, eps=1, min_samples=3, min_cluster_size=5)
                labels = sample.adata.obs[column]
                self.assertGreater(adjusted_rand_score([0]*20 + [1]*20, labels), .95)
                np.testing.assert_array_equal(sample.adata.X, before)
                self.assertEqual(sample.adata.uns['clustering']['algorithm'], algorithm)

    def test_density_noise_and_marker_edge_cases(self):
        sample = prepared_sample()
        sample.adata_processed.obsm['X_pca'][-1] = [100, 100, 100]
        column = sample.run_clustering(algorithm='dbscan', n_pcs=3, eps=1, min_samples=3)
        self.assertEqual(sample.adata.obs[column].iloc[-1], 'Noise')
        sample.compute_marker_genes()
        self.assertNotIn('Noise', sample.get_marker_genes())
        self.assertEqual(len(sample.get_marker_genes()), 2)
        for algorithm, options in [('dbscan', {'eps': 1e-10}), ('kmeans', {'n_clusters': 1})]:
            with self.subTest(algorithm=algorithm):
                sample.run_clustering(algorithm=algorithm, n_pcs=3, **options)
                sample.compute_marker_genes()
                self.assertEqual(sample.get_marker_genes(), {})

    def test_cache_tracks_method_parameters_and_marker_freshness(self):
        sample = prepared_sample()
        import src.utils.clustering as module
        with patch.object(module, 'fit_clusters', wraps=module.fit_clusters) as fitting:
            sample.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=2)
            sample.compute_marker_genes()
            self.assertTrue(sample.get_marker_genes())
            sample.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=2, eps=4)
            self.assertEqual(fitting.call_count, 1)
            sample.run_clustering(algorithm='kmeans', n_pcs=3, n_clusters=3)
            self.assertEqual(fitting.call_count, 2)
            self.assertEqual(sample.get_marker_genes(), {})
            sample.run_clustering(algorithm='agglomerative', n_pcs=3, n_clusters=3)
            self.assertEqual(fitting.call_count, 3)

    def test_aggregate_algorithms_project_identical_labels_to_native_members(self):
        _, _, _, bins = fixture()
        aggregate = bins._get_aggregated_sample(2)
        aggregate.adata_processed = aggregate.adata.copy()
        aggregate.adata_processed.obsm['X_pca'] = np.arange(aggregate.adata.n_obs*2).reshape(-1, 2)
        aggregate._step_params['pca'] = 2
        column = bins.run_clustering(algorithm='kmeans', n_clusters=2, n_pcs=2,
                    resolution_level=bins.get_resolution_level_for_factor(2), compute_markers=False)
        projected = bins.aggregation_names(2).map(aggregate.adata.obs[aggregate.cluster_col].astype(str))
        self.assertEqual(bins.adata.obs[column].tolist(), projected.tolist())
        self.assertEqual(column, 'kmeans_2')
        self.assertEqual(bins._clustering_factor, 2)

    def test_ui_parameter_visibility_rebuild_and_serialization(self):
        layer = make_layer()
        manager = layer.clustering_manager
        doc = Document()
        manager.get_widgets(doc)
        widgets = manager.get_widgets(doc)
        controls = manager._algorithm_controls
        self.assertEqual(len(controls['algorithm']._callbacks['value']), 1)
        for algorithm, fields in {
            'leiden': set(), 'kmeans': {'n_clusters'}, 'dbscan': {'eps', 'min_samples'},
            'hdbscan': {'min_cluster_size', 'min_samples'},
            'gaussian_mixture': {'n_clusters', 'covariance_type'},
            'agglomerative': {'n_clusters', 'linkage'},
        }.items():
            controls['algorithm'].value = algorithm
            visible = {key for key, widget in controls.items() if key != 'algorithm' and widget.visible}
            self.assertEqual(visible, fields)
        from bokeh.layouts import column
        doc.add_root(column(*widgets))
        doc.to_json()
        self.assertEqual(len(controls['algorithm'].options), 6)

    def test_selected_algorithm_reaches_run_and_noise_visibility_and_hover(self):
        from bokeh.events import ButtonClick
        from bokeh.models import Button
        from hd_viewer.board.utils.loading import LoadingLogger
        from tests.test_clustering_markers import ImmediateDoc
        layer = make_layer()
        manager = layer.clustering_manager
        widgets = manager.get_widgets(ImmediateDoc())
        manager._algorithm_controls['algorithm'].value = 'hdbscan'
        manager._algorithm_controls['min_cluster_size'].value = 12
        def work_only(logger, title, work, **callbacks):
            work()
        with patch.object(LoadingLogger, 'run', work_only), patch.object(manager, 'compute') as compute:
            button = next(w for w in widgets if isinstance(w, Button))
            button._trigger_event(ButtonClick(button))
            self.assertEqual(compute.call_args.kwargs['algorithm'], 'hdbscan')
            self.assertEqual(compute.call_args.kwargs['min_cluster_size'], 12)
        manager.cluster_col = 'dbscan_1'
        layer.sample.adata.obs['dbscan_1'] = ['Noise'] + ['0'] * (layer.sample.adata.n_obs - 1)
        manager._build_cluster_state()
        self.assertIn('Noise', manager._cluster_labels)
        manager._active_clusters = {'Noise'}
        layer.manager.source.data = {
            'obs_name': list(layer.sample.adata.obs_names),
            'dbscan_1': list(layer.sample.adata.obs['dbscan_1']),
        }
        manager._apply_visibility_mask()
        self.assertEqual(sum(layer._view), 1)
        layer.manager.source.data = {'obs_name': ['b0'], 'dbscan_1': ['Noise']}
        self.assertIn('Cluster: @{dbscan_1}', layer.get_hover_tooltips())

    def test_invalid_algorithm_and_cluster_counts(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported'):
            prepared_sample().run_clustering(algorithm='invalid')
        for count in (0, 50, 2.5):
            with self.assertRaisesRegex(ValueError, 'Number of clusters'):
                fit_clusters(np.ones((4, 3)), 'kmeans', n_clusters=count)


if __name__ == '__main__':
    unittest.main()
