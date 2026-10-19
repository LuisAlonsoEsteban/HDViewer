import threading
import time
import unittest
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from bokeh.models import Div

from bokeh_server.components.utils.groups import GROUPS_COL, GroupsManager


class ImmediateDoc:
    def add_next_tick_callback(self, callback):
        callback()


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
                                alias="Spot", _figure_de_div=Div(text=""))
    manager = GroupsManager(component)
    manager._clusters = {"A": [("g0", ">", "5", "AND")], "B": [("g1", ">", "5", "AND")]}
    manager._de_selected = {"A", "B"}
    # Only A is shown on the map, so the map's group column does not contain B
    manager._active_clusters = {"A"}
    adata.obs[GROUPS_COL] = np.where(adata.obs["grp_A"], "A", None)
    return manager, component


class GroupsDETests(unittest.TestCase):
    def test_between_groups_de_uses_hidden_groups_too(self):
        manager, component = make_manager()
        manager._de_status_div = Div(text="")
        done = threading.Event()
        component._figure_de_div.on_change("text", lambda attr, old, new: done.set())

        manager._run_de_thread(["A", "B"], mode=1, doc=ImmediateDoc())
        deadline = time.time() + 120
        while not done.is_set() and time.time() < deadline:
            time.sleep(0.2)

        self.assertTrue(done.is_set(), manager._de_status_div.text)
        self.assertEqual(manager._de_header, "DE: A vs B")
        self.assertIn("g0", component._figure_de_div.text)
        self.assertIn("g1", component._figure_de_div.text)

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
        component._figure_de_div.text = "results"

        manager._active_clusters = set()             # hidden on the map
        manager._clear_de_if_cluster_deleted()
        self.assertEqual(component._figure_de_div.text, "results")

        del manager._clusters["B"]
        manager._clear_de_if_cluster_deleted()
        self.assertEqual(component._figure_de_div.text, "")


if __name__ == "__main__":
    unittest.main()
