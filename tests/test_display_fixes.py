"""Cluster label colours and embedding colours follow the layer's assigned colours."""
import unittest

from hd_viewer.core import embedding
from tests.test_display_controls import clustered_layer


class ClusterLabelTests(unittest.TestCase):
    def test_each_cluster_is_listed_in_its_assigned_color(self):
        view, _ = clustered_layer()
        clustering = view.clustering
        summary = clustering.cluster_summary()
        self.assertEqual([entry["label"] for entry in summary], clustering._cluster_labels)
        for entry in summary:
            self.assertEqual(entry["colour"], clustering._color_map[entry["label"]])


class EmbeddingColorTests(unittest.TestCase):
    def setUp(self):
        self.view, self.cluster_col = clustered_layer()
        self.sample = self.view.sample
        self.names = list(self.sample.adata.obs_names.astype(str))

    def test_clusters_use_their_assigned_colors_for_every_observation(self):
        cmap = self.view.clustering._color_map
        colors, labels = embedding.colours(self.sample, self.cluster_col, self.names, palette=cmap)
        self.assertEqual(len(colors), len(self.names))
        self.assertEqual(colors, [cmap[label] for label in labels])

    def test_genes_use_a_gradient_on_the_layer_scale(self):
        self.sample.expression_scale = "log1p"
        colors, values = embedding.colours(self.sample, "g0", self.names)
        self.assertGreater(len(set(colors)), 2)
        self.assertNotIn(embedding.DEFAULT_COLOUR, colors)
        raw, _ = embedding.colours(self.sample, "", self.names)
        self.assertEqual(set(raw), {embedding.DEFAULT_COLOUR})


if __name__ == "__main__":
    unittest.main()
