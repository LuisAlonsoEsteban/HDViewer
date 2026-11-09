"""Clustering and embeddings of a sample share its PCA and neighbor graph."""
import unittest
from unittest.mock import patch

import scanpy as sc

from hd_viewer.core import embedding
from hd_viewer.settings.config import embedding_kwargs
from tests.test_clustering_markers import make_layer


class SharedStepTests(unittest.TestCase):
    def test_each_step_runs_once_across_clustering_and_embeddings(self):
        calls = {"pca": 0, "neighbors": 0, "umap": 0}

        def counted(name, function):
            def wrapper(*args, **kwargs):
                calls[name] += 1
                return function(*args, **kwargs)
            return wrapper

        with patch.object(sc.pp, "pca", counted("pca", sc.pp.pca)), \
                patch.object(sc.pp, "neighbors", counted("neighbors", sc.pp.neighbors)), \
                patch.object(sc.tl, "umap", counted("umap", sc.tl.umap)):
            view = make_layer()
            clustering, sample = view.clustering, view.sample
            options = embedding_kwargs()              # as the embedding API runs them
            clustering.compute(1.0)
            embedding.frame(sample, clustering, "pca", **options)    # reuses the clustering's PCA
            embedding.frame(sample, clustering, "umap", **options)   # reuses its neighbor graph
            clustering.compute(0.5)                                  # only Leiden is recomputed
            embedding.frame(sample, clustering, "umap", **options)   # still cached
        self.assertEqual(calls, {"pca": 1, "neighbors": 1, "umap": 1})

    def test_markers_of_replaced_labels_are_not_shown(self):
        clustering = make_layer().clustering
        clustering.compute(1.0)
        clustering.compute_markers()
        self.assertTrue(clustering.get_marker_genes())
        clustering.compute(2.0)
        self.assertEqual(clustering.get_marker_genes(), {})
        clustering.compute_markers()
        self.assertTrue(clustering.get_marker_genes())


if __name__ == "__main__":
    unittest.main()
