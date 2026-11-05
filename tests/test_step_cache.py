"""Clustering and embeddings of a sample share its PCA and neighbor graph."""
import unittest
from unittest.mock import patch

import scanpy as sc

from hd_viewer.board.layout.embedding import EmbeddingPanel
from tests.test_clustering_markers import ImmediateDoc, make_layer


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
            layer = make_layer()
            clustering = layer.clustering_manager
            panel = EmbeddingPanel(ImmediateDoc())
            panel.set_components([layer])
            clustering.compute(1.0)
            panel.compute_embedding(layer, "pca")          # reuses the clustering's PCA
            panel.compute_embedding(layer, "umap")         # reuses its neighbor graph
            clustering.compute(0.5)                        # only Leiden is recomputed
            panel.compute_embedding(layer, "umap")         # still cached
        self.assertEqual(calls, {"pca": 1, "neighbors": 1, "umap": 1})

    def test_markers_of_replaced_labels_are_not_shown(self):
        clustering = make_layer().clustering_manager
        clustering.compute(1.0)
        clustering.compute_markers()
        self.assertTrue(clustering.get_marker_genes())
        clustering.compute(2.0)
        self.assertEqual(clustering.get_marker_genes(), {})
        clustering.compute_markers()
        self.assertTrue(clustering.get_marker_genes())


if __name__ == "__main__":
    unittest.main()
