import threading
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.plotting import figure
from shapely.geometry import box

from bokeh_server.components.glyphs import Squares
from src.sample import Bins


class ImmediateDoc:
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


def make_bins(n=8):
    """n x n native bins of 1 x 1; gene g0 is high in the left half (x < n/2)."""
    names = [f"b{x}_{y}" for y in range(n) for x in range(n)]
    xs = np.array([x for y in range(n) for x in range(n)])
    index = pd.Index(names, name="barcode")
    gdf = gpd.GeoDataFrame({"index": names},
                           geometry=[box(x, y, x + 1, y + 1) for y in range(n) for x in range(n)],
                           index=index)
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(n * n, 20)).astype(float)
    X[xs < n // 2, 0] += 20
    adata = ad.AnnData(X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=[f"g{i}" for i in range(20)]))
    adata.obs["grp_A"] = xs < n // 2
    component = Squares(figure(), Bins(gdf, adata, pixels_per_bin=1, bin_size_um=2),
                        lazy_loading=True, alias="Bin")
    component._figure_de_layout = column()
    manager = component.groups_manager
    manager._clusters = {"A": []}
    manager._de_selected = {"A"}
    return component, manager, xs


class DELevelTests(unittest.TestCase):
    def test_membership_at_a_level_is_the_majority_of_its_bins(self):
        component, manager, xs = make_bins()
        cells = component.sample._get_aggregated_sample(2).adata.obs_names
        names = component.sample._aggregated_cell_names[2]

        # Cell of bins (0,0),(1,0),(0,1),(1,1): put 3 of its 4 bins in a set, then 1, then 2
        corner = names.index[names == names["b0_0"]]
        for k, expected in ((3, True), (1, False), (2, False)):
            mask = pd.Series(False, index=names.index)
            mask[corner[:k]] = True
            level = manager.de_analysis.masks_at_level({"S": mask.to_numpy()}, 2, cells)["S"]
            self.assertEqual(bool(level[list(cells).index(names["b0_0"])]), expected, k)

        # Left-half group: exactly the 8 cells of the left half
        level = manager.de_analysis.masks_at_level({"A": component.sample.adata.obs["grp_A"].to_numpy()}, 2, cells)["A"]
        self.assertEqual(int(level.sum()), 8)

    def test_de_runs_on_the_aggregated_cells(self):
        component, manager, _ = make_bins()
        done = threading.Event()
        component._figure_de_layout.on_change("children", lambda attr, old, new: done.set())
        manager._run_de_thread(["A"], 0, ImmediateDoc(), factor=2)
        self.assertTrue(done.wait(120), manager._de_status_div.text)

        self.assertEqual(manager._de_header, "DE vs Rest · 2 × 2 bins")
        self.assertEqual(manager.de_analysis.prepared_data[2].n_obs, 16)      # 2x2 cells, not 64 bins
        top = manager._de_results["A"].sort_values("scores", ascending=False)["names"].iloc[0]
        self.assertEqual(top, "g0")

    def test_resolution_choice_defaults_to_the_displayed_level(self):
        component, manager, _ = make_bins()
        widgets = manager.get_de_widgets(ImmediateDoc())
        select = next(w for w in widgets if getattr(w, "title", None) == "Resolution")
        self.assertEqual(select.value, "displayed")
        self.assertEqual(manager._de_factor_for(select), 1)              # zoom-driven, native

        component.manager.set_aggregation_factor(4)                     # bin panel: fixed 4 x 4
        self.assertEqual(manager._de_factor_for(select), 4)
        self.assertEqual(select.options[0], ("displayed", "Displayed level (4 × 4)"))
        select.value = "2"
        self.assertEqual(manager._de_factor_for(select), 2)

    def test_spots_have_no_resolution_choice(self):
        from tests.test_groups_de import make_manager
        manager, _ = make_manager()
        self.assertIsNone(manager._de_level_select())


if __name__ == "__main__":
    unittest.main()
