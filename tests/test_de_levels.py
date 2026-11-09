"""Differential expression at a bin aggregation level."""
import time
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.analysis.de import DEAnalysis, DEPresenter
from hd_viewer.startup.sample import Bins
from tests.test_clustering_markers import api_client
from tests.test_groups_de import make_sample as make_spots


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
    return Bins(gdf, adata, pixels_per_bin=1, bin_size=2), xs


class DELevelTests(unittest.TestCase):
    def test_membership_at_a_level_is_the_majority_of_its_bins(self):
        sample, xs = make_bins()
        analysis = DEAnalysis(sample)
        cells = sample._get_aggregated_sample(2).adata.obs_names
        names = sample._aggregated_cell_names[2]

        # Cell of bins (0,0),(1,0),(0,1),(1,1): put 3 of its 4 bins in a set, then 1, then 2
        corner = names.index[names == names["b0_0"]]
        self.assertEqual(len(corner), 4)
        for k, expected in ((3, True), (1, False), (2, False)):
            mask = pd.Series(False, index=names.index)
            mask[corner[:k]] = True
            level = analysis.masks_at_level({"S": mask.to_numpy()}, 2, cells)["S"]
            self.assertEqual(bool(level[list(cells).index(names["b0_0"])]), expected, k)

        # Left-half group: exactly the 8 cells of the left half
        level = analysis.masks_at_level({"A": xs < 4}, 2, cells)["A"]
        self.assertEqual(int(level.sum()), 8)

    def test_de_runs_on_the_aggregated_cells(self):
        sample, xs = make_bins()
        analysis = DEAnalysis(sample)
        result = analysis.run({"A": xs < 4}, 0, source="groups", factor=2)
        prepared = DEPresenter(alias="Bin", sample_type="bins").prepare(result)

        self.assertEqual(prepared.header, "DE vs Rest · 2 × 2 bins")
        self.assertEqual(analysis.prepared_data[2].n_obs, 16)      # 4x4 cells, not 64 bins
        self.assertEqual(result.factor, 2)
        top = result.results["A"].sort_values("scores", ascending=False)["names"].iloc[0]
        self.assertEqual(top, "g0")
        np.testing.assert_array_equal(sample.adata.X[xs < 4, 0] >= 20, True)   # counts untouched

    def test_resolution_choice_defaults_to_the_displayed_level(self):
        sample, xs = make_bins()
        client, session = api_client(self, sample)
        view = session.views["bins"]
        client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=-1, x1=9, y0=-1, y1=9, width=900, height=900))
        options = client.get("/api/layers/bins/de/options").json()
        self.assertEqual(options["displayed_factor"], 1)                 # zoom-driven, native
        self.assertEqual(options["factors"], list(sample.AGGREGATION_FACTORS))

        view.set_aggregation(4)                                          # bin panel: fixed 4 x 4
        self.assertEqual(client.get("/api/layers/bins/de/options").json()["displayed_factor"], 4)

        view.groups.add_fixed_group("A", xs < 4)

        def run(**body):
            self.assertEqual(client.post("/api/layers/bins/de", json=dict(groups=["A"], **body)).status_code, 202)
            for _ in range(2400):
                if view.job["state"] != "running":
                    return view.job
                time.sleep(.05)
            self.fail("DE did not finish")

        job = run()                               # the displayed level: 2 x 2 cells of 4 x 4 bins
        self.assertEqual(job["state"], "failed")
        self.assertIn("≥3 4 × 4 bins", job["message"])
        self.assertEqual(run(factor=2)["state"], "done")                  # an explicit choice
        self.assertEqual(view.de["header"], "DE vs Rest · 2 × 2 bins")

    def test_spots_have_no_resolution_choice(self):
        client, _ = api_client(self, make_spots(), kind="spots", layer_id="spots")
        options = client.get("/api/layers/spots/de/options").json()
        self.assertEqual(options["factors"], [])
        self.assertEqual(options["displayed_factor"], 1)


if __name__ == "__main__":
    unittest.main()
