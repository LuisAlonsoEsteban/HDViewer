import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer import settings
from hd_viewer.core.groups import GROUPS_COL
from hd_viewer.core.palette import CLUSTER_PALETTE, GROUP_PALETTE, palette_color
from hd_viewer.startup.sample import Bins
from tests.fixtures import make_view


def make_view_with_clusters(n=40):
    index = pd.Index([f"b{i}" for i in range(n)], name="barcode")
    gdf = gpd.GeoDataFrame({"index": list(index)}, geometry=[box(i, 0, i + 1, 1) for i in range(n)],
                           index=index)
    adata = ad.AnnData(np.ones((n, 2)), obs=pd.DataFrame(index=index.astype(str)),
                       var=pd.DataFrame(index=["g0", "g1"]))
    adata.obs["leiden"] = pd.Categorical([str(i % 25) for i in range(n)])
    adata.obs["some_metadata"] = pd.Categorical(["x", "y"] * (n // 2))
    view = make_view(Bins(gdf, adata, pixels_per_bin=1, bin_size=16), kind="bins")
    view.update_viewport(-1, n + 1, -1, 2, 900, 900)
    return view


class PaletteTests(unittest.TestCase):
    def test_groups_and_clusters_never_share_a_colour(self):
        self.assertFalse({c.lower() for c in CLUSTER_PALETTE} & {c.lower() for c in GROUP_PALETTE})
        self.assertNotEqual(palette_color(GROUP_PALETTE, 0), palette_color(CLUSTER_PALETTE, 0))

    def test_more_than_20_clusters_get_distinct_colours(self):
        colours = [palette_color(CLUSTER_PALETTE, i) for i in range(len(CLUSTER_PALETTE))]
        self.assertGreater(len(colours), 50)
        self.assertEqual(len(set(colours)), len(colours))

    def test_the_map_uses_the_assigned_cluster_colours(self):
        view = make_view_with_clusters()
        clustering = view.clustering
        clustering.cluster_col = "leiden"
        clustering._build_cluster_state()
        view.colour_by_categories("leiden", clustering._color_map, "clusters")
        param = view.rows.get_param("color")

        self.assertEqual(param.category_colors, clustering._color_map)
        # Category values come in any order (here 7, 0, 24): colours follow the label
        colours = param.decode(param.get_values(pd.Series(["7", "0", "24"])))
        self.assertEqual(colours, [clustering._color_map[c] for c in ("7", "0", "24")])
        self.assertEqual(len(set(clustering._color_map.values())), 25)
        self.assertEqual(clustering._color_map["0"], CLUSTER_PALETTE[0])
        # The loaded rows too
        labels = view.sample.adata.obs["leiden"].astype(str)
        shown = param.decode(view.rows.sink.data["color"])
        names = list(view.rows.loaded_obs_names)
        self.assertEqual(shown, [clustering._color_map[labels[name]] for name in names])

    def test_the_map_uses_the_assigned_group_colours(self):
        view = make_view_with_clusters()
        groups = view.groups
        name = groups.add_fixed_group("Tumour", np.arange(40) < 10)     # as created in the Groups tab
        groups.refresh_active_labels()
        self.assertEqual(groups._color_map[name], GROUP_PALETTE[0])
        view.colour_by_categories(GROUPS_COL, groups._color_map, "groups")
        self.assertEqual(view.rows.get_param("color").category_colors["Tumour"], GROUP_PALETTE[0])

    def test_other_categorical_columns_keep_automatic_colours(self):
        view = make_view_with_clusters()
        view.colour_by_column("some_metadata")
        param = view.rows.get_param("color")
        self.assertEqual(view.colour_by, dict(kind="category", name="some_metadata"))
        # No assigned colours: the automatic categorical palette, in order of appearance
        self.assertEqual(param.category_colors, {"x": settings.CATEGORICAL_PALETTE[0],
                                                 "y": settings.CATEGORICAL_PALETTE[1]})


if __name__ == "__main__":
    unittest.main()
