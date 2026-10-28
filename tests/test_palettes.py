import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from shapely.geometry import box

from bokeh_server.components.glyphs import Squares
from bokeh_server.components.utils.groups import GROUPS_COL
from bokeh_server.components.utils.palettes import CLUSTER_PALETTE, GROUP_PALETTE, palette_color
from src.sample import Bins


def make_component(n=40):
    index = pd.Index([f"b{i}" for i in range(n)], name="barcode")
    gdf = gpd.GeoDataFrame({"index": list(index)}, geometry=[box(i, 0, i + 1, 1) for i in range(n)],
                           index=index)
    adata = ad.AnnData(np.ones((n, 2)), obs=pd.DataFrame(index=index.astype(str)),
                       var=pd.DataFrame(index=["g0", "g1"]))
    adata.obs["leiden"] = pd.Categorical([str(i % 25) for i in range(n)])
    return Squares(figure(), Bins(gdf, adata, pixels_per_bin=1, bin_size=16),
                   lazy_loading=True, alias="Bin")


class PaletteTests(unittest.TestCase):
    def test_groups_and_clusters_never_share_a_colour(self):
        self.assertFalse({c.lower() for c in CLUSTER_PALETTE} & {c.lower() for c in GROUP_PALETTE})
        self.assertNotEqual(palette_color(GROUP_PALETTE, 0), palette_color(CLUSTER_PALETTE, 0))

    def test_more_than_20_clusters_get_distinct_colours(self):
        colours = [palette_color(CLUSTER_PALETTE, i) for i in range(len(CLUSTER_PALETTE))]
        self.assertGreater(len(colours), 50)
        self.assertEqual(len(set(colours)), len(colours))

    def test_the_map_uses_the_assigned_cluster_colours(self):
        component = make_component()
        clustering = component.clustering_manager
        clustering.cluster_col = "leiden"
        clustering._build_cluster_state()
        component.add_parameter(attr="color", colname="leiden", static=False, not_sequenced_val="#ccc")
        param = component.manager.get_param("color")

        self.assertEqual(param.category_colors, clustering._color_map)
        # Category values come in any order (here 7, 0, 24): colours follow the label
        colours = param.get_values(pd.Series(["7", "0", "24"]))
        self.assertEqual(colours, [clustering._color_map[c] for c in ("7", "0", "24")])
        self.assertEqual(len(set(clustering._color_map.values())), 25)

    def test_the_map_uses_the_assigned_group_colours(self):
        component = make_component()
        component.sample.adata.obs[GROUPS_COL] = "Tumour"   # as created in the Groups tab
        component.groups_manager._color_map = {"Tumour": palette_color(GROUP_PALETTE, 0)}
        component.add_parameter(attr="color", colname=GROUPS_COL, static=False, not_sequenced_val="#ccc")
        self.assertEqual(component.manager.get_param("color").category_colors["Tumour"],
                         GROUP_PALETTE[0])

    def test_other_categorical_columns_keep_automatic_colours(self):
        component = make_component()
        self.assertEqual(component.known_category_colors("some_metadata"), {})


if __name__ == "__main__":
    unittest.main()
