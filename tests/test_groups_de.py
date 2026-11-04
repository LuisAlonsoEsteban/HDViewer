import threading
import time
import unittest
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.models import Div
from bokeh.layouts import column
from bokeh.document import Document
from bokeh.models import HoverTool, Plot

from hd_viewer.board.layout.layer.groups import GROUPS_COL, GroupsManager


class ImmediateDoc:
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


def make_manager():
    """60 spots: group A = 0-19 (gene g0 high), group B = 20-39 (gene g1 high)."""
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(60, 30)).astype(float)
    X[:20, 0] += 30
    X[20:40, 1] += 30
    adata = ad.AnnData(X, obs=pd.DataFrame(index=[f"s{i}" for i in range(60)]),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(30)]))
    adata.obs["grp_A"] = np.arange(60) < 20
    adata.obs["grp_B"] = (np.arange(60) >= 20) & (np.arange(60) < 40)

    component = SimpleNamespace(manager=SimpleNamespace(sample=SimpleNamespace(adata=adata)),
                                alias="Spot", _figure_de_layout=column())
    manager = GroupsManager(component)
    manager._clusters = {"A": [("g0", ">", "5", "AND")], "B": [("g1", ">", "5", "AND")]}
    manager._de_selected = {"A", "B"}
    # Only A is shown on the map, so the map's group column does not contain B
    manager._active_clusters = {"A"}
    adata.obs[GROUPS_COL] = np.where(adata.obs["grp_A"], "A", None)
    return manager, component


def result_html(component):
    return "".join(model.text for model in component._figure_de_layout.references() if isinstance(model, Div))


class GroupsDETests(unittest.TestCase):
    def test_results_attach_live_figures_and_clear_them_on_group_deletion(self):
        manager, component = make_manager()
        document = Document()
        document.add_root(component._figure_de_layout)
        self.assertTrue(run_and_wait(manager, component, ["A", "B"], 1))
        plots = list(component._figure_de_layout.select({'type': Plot}))
        self.assertEqual(len(plots), 2)
        self.assertTrue(all(plot.document is document for plot in plots))
        self.assertTrue(all(any(isinstance(tool, HoverTool) for tool in plot.tools) for plot in plots))
        document.to_json()
        del manager._clusters['B']
        manager._clear_de_if_cluster_deleted()
        self.assertEqual(component._figure_de_layout.children, [])
        self.assertFalse(any(plot.id in document.models for plot in plots))

    def test_between_groups_de_uses_hidden_groups_too(self):
        manager, component = make_manager()
        manager._de_status_div = Div(text="")
        done = threading.Event()
        component._figure_de_layout.on_change("children", lambda attr, old, new: done.set())

        manager._run_de_thread(["A", "B"], mode=1, doc=ImmediateDoc())
        deadline = time.time() + 120
        while not done.is_set() and time.time() < deadline:
            time.sleep(0.2)

        self.assertTrue(done.is_set(), manager._de_status_div.text)
        self.assertEqual(manager._de_header, "DE: A vs B")
        self.assertIn("g0", result_html(component))
        self.assertIn("g1", result_html(component))

    def test_de_group_list_follows_created_and_deleted_groups(self):
        manager, _ = make_manager()
        widgets = manager.get_de_widgets(ImmediateDoc())
        checkbox = manager._de_groups_checkbox
        self.assertIn(checkbox, widgets)
        self.assertEqual(checkbox.labels, ["A", "B"])
        self.assertEqual(checkbox.active, [0, 1])

        checkbox.active = [1]                       # compare only B
        self.assertEqual(manager._de_selected, {"B"})
        manager._clusters["C"] = []
        manager._de_selected.add("C")                # new groups start selected
        manager._refresh_de_groups()
        self.assertEqual(checkbox.labels, ["A", "B", "C"])
        self.assertEqual(checkbox.active, [1, 2])

        del manager._clusters["B"]
        manager._refresh_de_groups()
        self.assertEqual(checkbox.labels, ["A", "C"])
        self.assertEqual(checkbox.active, [1])

    def test_hiding_a_group_keeps_de_results_but_deleting_it_clears_them(self):
        manager, component = make_manager()
        manager._de_status_div = Div(text="")
        manager._de_clusters = ["A", "B"]
        component._figure_de_layout.children = [Div(text="results")]

        manager._active_clusters = set()             # hidden on the map
        manager._clear_de_if_cluster_deleted()
        self.assertEqual(result_html(component), "results")

        del manager._clusters["B"]
        manager._clear_de_if_cluster_deleted()
        self.assertEqual(result_html(component), "")


def with_clusters(manager):
    """Leiden-like clusters: '0' = spots 0-29 (g0 high), '1' = 30-59 (g2 high)."""
    adata = manager.sample.adata
    adata.obs["leiden"] = np.where(np.arange(60) < 30, "0", "1")
    adata.X[30:, 2] += 30
    manager.component.clustering_manager = SimpleNamespace(
        is_clustered=True, cluster_col="leiden", _cluster_labels=["0", "1"],
        _color_map={"0": "#123456", "1": "#abcdef"})
    return manager


def run_and_wait(manager, component, *args, **kwargs):
    manager._de_status_div = manager._de_status_div or Div(text="")
    done = threading.Event()
    component._figure_de_layout.on_change("children", lambda attr, old, new: done.set())
    manager._run_de_thread(*args, doc=ImmediateDoc(), **kwargs)
    deadline = time.time() + 120
    while not done.is_set() and time.time() < deadline:
        time.sleep(0.2)
    return done.is_set()


class ClusterDETests(unittest.TestCase):
    def test_groups_and_clusters_are_listed_side_by_side(self):
        manager, _ = make_manager()
        widgets = manager.get_de_widgets(ImmediateDoc())
        mode = next(w for w in widgets if getattr(w, "labels", None) == ["vs Rest", "Between selected"])
        self.assertNotIn("Entre clusters", mode.labels)
        self.assertEqual(manager._de_groups_checkbox.labels, ["A", "B"])
        self.assertEqual(manager._de_clusters_checkbox.labels, [])   # no clustering yet
        self.assertIn("Run clustering in the Clustering tab", manager._de_clusters_hint.text)

        with_clusters(manager).refresh_de_clusters()
        self.assertEqual(manager._de_clusters_checkbox.labels, ["Cluster 0", "Cluster 1"])
        # Unticked, so that one group can be compared with one cluster
        self.assertEqual(manager._de_clusters_checkbox.active, [])
        self.assertEqual(manager._de_groups_checkbox.active, [0, 1])

    def test_a_group_named_like_a_cluster_stays_distinct(self):
        manager, _ = make_manager()
        manager._clusters["Cluster 0"] = []
        with_clusters(manager)
        self.assertEqual(manager._de_options("clusters"), ["Cluster 0 (clustering)", "Cluster 1"])
        mask = manager._de_masks("clusters", ["Cluster 0 (clustering)"])["Cluster 0 (clustering)"]
        self.assertEqual(int(mask.sum()), 30)

    def test_group_vs_cluster_leaves_out_their_overlap(self):
        manager, component = make_manager()
        with_clusters(manager)
        # Group B = spots 20-39, Cluster 1 = spots 30-59: spots 30-39 are in both
        masks = {**manager._de_masks("groups", ["B"]), **manager._de_masks("clusters", ["Cluster 1"])}
        self.assertTrue(run_and_wait(manager, component, ["B", "Cluster 1"], 1,
                                     masks=masks, source="mixed"))
        self.assertEqual(manager._de_header, "DE: B vs Cluster 1")
        top = {name: df.sort_values("scores", ascending=False)["names"].iloc[0]
               for name, df in manager._de_results.items()}
        # B outside the overlap (20-29) has g1 high; Cluster 1 outside it (40-59) has g2 high
        self.assertEqual(top, {"B": "g1", "Cluster 1": "g2"})
        self.assertIn("10 overlapping spots left out", manager._de_status_div.text)

    def test_run_button_compares_a_ticked_group_with_a_ticked_cluster(self):
        from bokeh.events import ButtonClick
        manager, component = make_manager()
        with_clusters(manager)
        widgets = manager.get_de_widgets(ImmediateDoc())
        mode = next(w for w in widgets if getattr(w, "labels", None) == ["vs Rest", "Between selected"])
        button = next(w for w in widgets if getattr(w, "label", None) == "▶ Run DE")
        manager._de_groups_checkbox.active = [1]          # B
        manager._de_clusters_checkbox.active = [1]        # Cluster 1
        mode.active = 1

        done = threading.Event()
        component._figure_de_layout.on_change("children", lambda attr, old, new: done.set())
        button._trigger_event(ButtonClick(button))
        self.assertTrue(done.wait(120), manager._de_status_div.text)
        self.assertEqual(set(manager._de_results), {"B", "Cluster 1"})
        self.assertEqual(manager._de_results_source, "mixed")
        self.assertIn("#abcdef", result_html(component))   # Cluster 1's colour

    def test_between_clusters_de(self):
        manager, component = make_manager()
        with_clusters(manager)
        names = ["Cluster 0", "Cluster 1"]
        self.assertTrue(run_and_wait(
            manager, component, names, 1, masks=manager._de_masks("clusters", names),
            colors=manager._de_colors("clusters"), source="clusters"))
        self.assertEqual(manager._de_header, "DE: Cluster 0 vs Cluster 1")
        top = {name: df.sort_values("scores", ascending=False)["names"].iloc[0]
               for name, df in manager._de_results.items()}
        self.assertEqual(top, {"Cluster 0": "g0", "Cluster 1": "g2"})
        self.assertIn("#123456", result_html(component))   # cluster colours
        # Same status line as the other long jobs: finished, with the elapsed time
        status = manager._de_status_div.text
        self.assertIn("DE done: 2 clusters, results", status)
        self.assertRegex(status, r"\(\d+\.\ds\)")

    def test_failure_is_reported_in_the_status_line(self):
        manager, component = make_manager()
        manager._de_logger = None
        logger = manager._get_de_logger(ImmediateDoc())
        finished = threading.Event()
        logger.div.on_change("text", lambda attr, old, new: "Error" in new and finished.set())
        # Group A has ≥3 spots but nothing else is selected: no comparison possible
        manager._run_de_thread(["A"], 1, ImmediateDoc())
        self.assertTrue(finished.wait(60), logger.div.text)
        self.assertIn("needs at least 2 groups", logger.div.text)
        self.assertFalse(logger.running)

    def test_reclustering_clears_cluster_results_only(self):
        manager, component = make_manager()
        with_clusters(manager)
        manager._de_status_div = Div(text="")
        manager._de_clusters, manager._de_results_source = ["A"], "groups"
        component._figure_de_layout.children = [Div(text="group results")]
        manager.refresh_de_clusters()
        self.assertEqual(result_html(component), "group results")

        manager._de_clusters, manager._de_results_source = ["Cluster 0"], "clusters"
        manager.refresh_de_clusters()
        self.assertEqual(result_html(component), "")


if __name__ == "__main__":
    unittest.main()
