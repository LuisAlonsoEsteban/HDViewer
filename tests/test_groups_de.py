"""Differential expression between custom groups and clusters (core sets + HTTP API).

60 spots: group A = spots 0-19 (gene g0 high), group B = 20-39 (gene g1 high);
with clusters, '0' = spots 0-29 and '1' = spots 30-59 (g2 high).
"""
import json
import time
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.analysis.de import DEAnalysis, DEPresenter
from hd_viewer.core.de_sets import de_colours, de_inputs, de_masks, de_options
from hd_viewer.startup.sample import Spots
from tests.fixtures import make_view
from tests.test_clustering_markers import api_client

CONDITIONS = {"A": [dict(column="g0", op=">", value="5")], "B": [dict(column="g1", op=">", value="5")]}


def make_sample(clusters=False):
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(60, 30)).astype(float)
    X[:20, 0] += 30
    X[20:40, 1] += 30
    names = [f"s{i}" for i in range(60)]
    if clusters:
        X[30:, 2] += 30
    gdf = gpd.GeoDataFrame({"index": names}, geometry=[box(i, 0, i + 1, 1) for i in range(60)],
                           index=pd.Index(names, name="barcode"))
    adata = ad.AnnData(X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=[f"g{i}" for i in range(30)]))
    if clusters:
        adata.obs["leiden"] = np.where(np.arange(60) < 30, "0", "1")
    return Spots(gdf, adata)


def with_clusters(view):
    """Leiden-like clusters (on a `make_sample(clusters=True)` sample), as a finished clustering job."""
    clustering = view.clustering
    clustering._cluster_columns["leiden"] = clustering.cluster_col = "leiden"
    clustering._build_cluster_state()
    clustering._color_map.update({"0": "#123456", "1": "#abcdef"})
    return view


def make_groups_view(clusters=False):
    """View with groups A and B; only A is shown on the map."""
    view = make_view(make_sample(clusters), kind="spots")
    groups = view.groups
    for name, conditions in CONDITIONS.items():
        groups.commit_condition_group(*groups.prepare_condition_group(
            name, [(c["column"], c["op"], c["value"], "AND") for c in conditions]))
    groups._active_groups.discard("B")
    groups.refresh_active_labels()
    return with_clusters(view) if clusters else view


def run(view, chosen, mode):
    names, masks, colours, source = de_inputs(view.sample, chosen, view.groups, view.clustering)
    result = DEAnalysis(view.sample).run(masks, mode, source=source)
    return result, DEPresenter(alias="Spot", colors=colours).prepare(result)


def top_genes(result):
    return {name: df.sort_values("scores", ascending=False)["names"].iloc[0] for name, df in result.results.items()}


class Api:
    """The 60-spot layer behind the HTTP API, with groups A and B (B hidden)."""
    def __init__(self, test, clusters=False):
        self.test = test
        self.client, self.session = api_client(test, make_sample(clusters), kind="spots", layer_id="spots")
        self.view = self.session.views["spots"]
        for name, conditions in CONDITIONS.items():
            response = self.client.post("/api/layers/spots/groups", json=dict(name=name, conditions=conditions))
            test.assertEqual(response.status_code, 200, response.text)
        test.assertEqual(self.client.patch("/api/layers/spots/groups/B", json=dict(visible=False)).status_code, 200)
        if clusters:
            with_clusters(self.view)

    def get(self, path, **params):
        return self.client.get(f"/api/layers/spots/{path}", params=params)

    def run_de(self, **body):
        response = self.client.post("/api/layers/spots/de", json=body)
        self.test.assertEqual(response.status_code, 202, response.text)
        for _ in range(2400):
            job = self.view.job
            if job["state"] != "running":
                return dict(job)
            time.sleep(.05)
        self.test.fail("DE did not finish")


class GroupsDETests(unittest.TestCase):
    def test_results_have_live_plots_and_are_cleared_on_group_deletion_only(self):
        api = Api(self)
        job = api.run_de(groups=["A", "B"], mode=1)
        self.assertEqual(job["state"], "done", job["message"])
        plots = json.dumps(api.get("de/plots").json())
        self.assertIn('"doc"', plots)
        self.assertIn("HoverTool", plots)
        self.assertEqual(api.get("de/table").json()["sets"], ["A", "B"])

        # Hiding a group on the map keeps the results
        api.client.patch("/api/layers/spots/groups/A", json=dict(visible=False))
        self.assertTrue(api.get("de/options").json()["ready"])
        self.assertEqual(api.get("de/plots").status_code, 200)
        # Deleting a compared group clears them
        api.client.delete("/api/layers/spots/groups/B")
        self.assertFalse(api.get("de/options").json()["ready"])
        self.assertEqual(api.get("de/plots").status_code, 404)
        self.assertEqual(api.get("de/table").status_code, 404)

    def test_between_groups_de_uses_hidden_groups_too(self):
        view = make_groups_view()
        self.assertEqual(view.groups._active_groups, {"A"})       # B is not in the map's group column
        result, prepared = run(view, [("groups", "A"), ("groups", "B")], 1)
        self.assertEqual(prepared.header, "DE: A vs B")
        self.assertEqual(top_genes(result), {"A": "g0", "B": "g1"})
        self.assertIn("g0", prepared.table_html)
        self.assertIn("g1", prepared.table_html)

    def test_de_group_list_follows_created_and_deleted_groups(self):
        api = Api(self)
        names = lambda: [g["name"] for g in api.get("de/options").json()["groups"]]
        self.assertEqual(names(), ["A", "B"])
        api.client.post("/api/layers/spots/groups", json=dict(name="C", conditions=[
            dict(column="g2", op=">", value="0")]))
        self.assertEqual(names(), ["A", "B", "C"])
        api.client.delete("/api/layers/spots/groups/B")
        self.assertEqual(names(), ["A", "C"])
        colours = {g["name"]: g["colour"] for g in api.get("de/options").json()["groups"]}
        self.assertEqual(colours, {name: api.view.groups._color_map[name] for name in ("A", "C")})


class ClusterDETests(unittest.TestCase):
    def test_groups_and_clusters_are_listed_side_by_side(self):
        view = make_view(make_sample(clusters=True), kind="spots")
        groups = view.groups
        groups.commit_condition_group(*groups.prepare_condition_group("A", [("g0", ">", "5", "AND")]))
        self.assertEqual(de_options("groups", groups, view.clustering), ["A"])
        self.assertEqual(de_options("clusters", groups, view.clustering), [])   # no clustering yet
        with_clusters(view)
        self.assertEqual(de_options("clusters", groups, view.clustering), ["Cluster 0", "Cluster 1"])
        self.assertEqual(de_colours("clusters", groups, view.clustering),
                         {"Cluster 0": "#123456", "Cluster 1": "#abcdef"})

        from hd_viewer.api.routers.de import set_options
        options = set_options(view)
        self.assertEqual(options["clusters"], [dict(name="Cluster 0", colour="#123456"),
                                               dict(name="Cluster 1", colour="#abcdef")])
        self.assertEqual([g["name"] for g in options["groups"]], ["A"])

    def test_a_group_named_like_a_cluster_stays_distinct(self):
        view = make_groups_view(clusters=True)
        view.groups._groups["Cluster 0"] = []
        self.assertEqual(de_options("clusters", view.groups, view.clustering), ["Cluster 0 (clustering)", "Cluster 1"])
        mask = de_masks(view.sample, "clusters", ["Cluster 0 (clustering)"], view.groups,
                        view.clustering)["Cluster 0 (clustering)"]
        self.assertEqual(int(mask.sum()), 30)

    def test_group_vs_cluster_leaves_out_their_overlap(self):
        view = make_groups_view(clusters=True)
        # Group B = spots 20-39, Cluster 1 = spots 30-59: spots 30-39 are in both
        result, prepared = run(view, [("groups", "B"), ("clusters", "Cluster 1")], 1)
        self.assertEqual(result.source, "mixed")
        self.assertEqual(prepared.header, "DE: B vs Cluster 1")
        # B outside the overlap (20-29) has g1 high; Cluster 1 outside it (40-59) has g2 high
        self.assertEqual(top_genes(result), {"B": "g1", "Cluster 1": "g2"})
        self.assertEqual(result.n_overlap, 10)

    def test_api_compares_a_group_with_a_cluster(self):
        api = Api(self, clusters=True)
        job = api.run_de(groups=["B"], clusters=["Cluster 1"], mode=1)
        self.assertEqual(job["state"], "done", job["message"])
        self.assertIn("10 overlapping observations left out", job["message"])
        table = api.get("de/table").json()
        self.assertEqual(table["sets"], ["B", "Cluster 1"])
        self.assertEqual(table["colours"][1], "#abcdef")                # Cluster 1's colour
        self.assertEqual(table["colours"][0], api.view.groups._color_map["B"])
        self.assertEqual(api.view.de["compared"], ("mixed", ["B", "Cluster 1"]))

    def test_between_clusters_de(self):
        view = make_groups_view(clusters=True)
        names = ["Cluster 0", "Cluster 1"]
        result, prepared = run(view, [("clusters", name) for name in names], 1)
        self.assertEqual(prepared.header, "DE: Cluster 0 vs Cluster 1")
        self.assertEqual(top_genes(result), {"Cluster 0": "g0", "Cluster 1": "g2"})
        self.assertIn("#123456", prepared.table_html)   # cluster colours

        api = Api(self, clusters=True)
        job = api.run_de(clusters=names, mode=1)
        self.assertEqual(job["message"], "DE done: 2 clusters")
        self.assertEqual(api.get("de/table").json()["colours"], ["#123456", "#abcdef"])

    def test_failure_is_reported(self):
        view = make_groups_view()
        # Group A has ≥3 spots but nothing else is selected: no comparison possible
        masks = de_masks(view.sample, "groups", ["A"], view.groups)
        with self.assertRaisesRegex(ValueError, "needs at least 2 groups"):
            DEAnalysis(view.sample).run(masks, 1, source="groups")

        api = Api(self)
        response = api.client.post("/api/layers/spots/de", json=dict(groups=["A"], mode=1))
        self.assertEqual(response.status_code, 422)
        self.assertIn("needs at least 2", response.json()["detail"])
        # A failure in the job is reported in its status and leaves the layer free
        empty = np.zeros(60, dtype=bool)
        empty[:2] = True
        api.view.groups.add_fixed_group("Tiny", empty)
        job = api.run_de(groups=["Tiny", "A"], mode=1)
        self.assertEqual(job["state"], "failed")
        self.assertIn("needs at least 2 groups", job["message"])
        self.assertIsNone(api.view.de["running"])
        self.assertEqual(api.run_de(groups=["A"])["state"], "done")

    def test_reclustering_clears_cluster_results_only(self):
        view = make_groups_view(clusters=True)
        results = dict(header="DE vs Rest", merged=pd.DataFrame(), plot_data=[], names=["A"], colours={})
        view.de.update(results, compared=("groups", ["A"]))
        self.assertFalse(view.drop_comparisons(clusters=True))
        self.assertEqual(view.de["header"], "DE vs Rest")

        for source in ("clusters", "mixed"):
            with self.subTest(source=source):
                view.de.update(results, compared=(source, ["Cluster 0"]))
                self.assertTrue(view.drop_comparisons(clusters=True))
                self.assertIsNone(view.de["header"])
                self.assertIsNone(view.de["merged"])

        # A clustering finishing while a cluster DE runs discards that DE's result
        view.de.update(running=("clusters", ["Cluster 0"]), generation=3)
        view.drop_comparisons(clusters=True)
        self.assertEqual((view.de["running"], view.de["generation"]), (None, 4))


if __name__ == "__main__":
    unittest.main()
