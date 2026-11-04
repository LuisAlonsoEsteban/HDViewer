"""Toolbar stability, cluster label colors, readable menus and embedding colors."""
import re
import unittest
from pathlib import Path

import numpy as np
from bokeh.models import CheckboxGroup, HoverTool

from hd_viewer.board.layout.embedding import EmbeddingPanel
from tests.test_clustering_markers import ImmediateDoc, make_layer as make_clustered_layer
from tests.test_roi_manager import make_layer


def clustered_layer():
    layer = make_clustered_layer()
    cm = layer.clustering_manager
    cluster_col = cm.compute(1.0)
    cm.apply_result(cluster_col, ImmediateDoc())
    return layer, cluster_col


class ToolbarTests(unittest.TestCase):
    def test_parameter_changes_reuse_the_hover_tool(self):
        layer = make_layer("bins")
        hover, tools = layer.hover, list(layer.figure.tools)
        self.assertIsInstance(hover, HoverTool)
        layer.add_parameter(attr="line_width", colname="g1", colname_hovered=True)
        layer.remove_parameter("line_width")
        layer.refresh_hover()
        self.assertIs(layer.hover, hover)
        self.assertEqual(list(layer.figure.tools), tools)        # toolbar never rebuilt
        self.assertEqual(hover.renderers, [layer.glyph_renderer])


class ClusterLabelTests(unittest.TestCase):
    def test_each_cluster_is_listed_in_its_assigned_color(self):
        layer, _ = clustered_layer()
        cm = layer.clustering_manager
        boxes = [w for w in cm._build_cluster_control_widgets(ImmediateDoc()) if isinstance(w, CheckboxGroup)]
        self.assertEqual(len(boxes), len(cm._cluster_labels))
        for box, label in zip(boxes, cm._cluster_labels):
            css = " ".join(sheet.css for sheet in box.stylesheets if hasattr(sheet, "css"))
            self.assertIn(f"color: {cm._color_map[label]}", css)


class MenuStyleTests(unittest.TestCase):
    def test_autocomplete_suggestions_have_a_dark_background(self):
        css = (Path(__file__).resolve().parents[1] / "hd_viewer/static/css/widgets.css").read_text()
        rule = re.search(r"\.bk-menu \{([^}]*)\}", css)
        self.assertIsNotNone(rule)
        self.assertIn("background-color: #11231d", rule.group(1))


class EmbeddingColorTests(unittest.TestCase):
    def setUp(self):
        self.layer, self.cluster_col = clustered_layer()
        self.panel = EmbeddingPanel(ImmediateDoc())
        self.panel.set_components([self.layer])
        self.view = self.panel._component_views[self.layer.key]
        self.names = list(self.layer.sample.adata.obs_names.astype(str))

    def test_clusters_use_their_assigned_colors_for_every_observation(self):
        colors, labels = self.panel.colors_for(self.layer, self.cluster_col, self.names)
        cmap = self.layer.clustering_manager._color_map
        self.assertEqual(colors, [cmap[label] for label in labels])

    def test_genes_use_a_gradient_on_the_layer_scale(self):
        self.layer.sample.expression_scale = "log1p"
        colors, values = self.panel.colors_for(self.layer, "g0", self.names)
        self.assertGreater(len(set(colors)), 2)
        self.assertNotIn(self.panel.DEFAULT_COLOR, colors)
        raw, _ = self.panel.colors_for(self.layer, "", self.names)
        self.assertEqual(set(raw), {self.panel.DEFAULT_COLOR})

    def test_selector_recolors_and_follows_new_clusterings(self):
        self.view["source"].data = dict(obs_name=self.names, x=[0.0] * len(self.names),
                                        y=[0.0] * len(self.names), color=["red"] * len(self.names),
                                        value=[""] * len(self.names))
        self.view["color_input"].value = "g0"
        self.assertNotIn("red", self.view["source"].data["color"])
        cm = self.layer.clustering_manager
        cluster_col = cm.compute(2.0)
        cm.apply_result(cluster_col, ImmediateDoc())
        for callback in self.layer._clustering_listeners:
            callback(self.layer)
        self.assertEqual(self.view["color_input"].value, cluster_col)
        labels = self.view["source"].data["value"]
        self.assertEqual(self.view["source"].data["color"], [cm._color_map[label] for label in labels])
        self.assertIn(cluster_col, self.view["color_input"].completions)


if __name__ == "__main__":
    unittest.main()
