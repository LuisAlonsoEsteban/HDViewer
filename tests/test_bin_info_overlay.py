"""Bin resolution state shown by the bin-info overlay (drawn by the React front
end from `LayerView.bin_info()`, sent in every frame header): the displayed
aggregation factor and bin geometry, fixed vs zoom-driven levels and the
cluster minimum."""
import math
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Polygon

from hd_viewer.startup.sample import Bins
from tests.fixtures import loaded, make_view, show


class BinResolutionTests(unittest.TestCase):
    def setUp(self):
        # A native bin transformed by rotation followed by unequal axis scaling.
        angle = math.radians(30)
        rotation = np.array([[math.cos(angle), -math.sin(angle)],
                             [math.sin(angle), math.cos(angle)]])
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) @ rotation.T
        corners = corners * [1, 2] + [10, 10]
        self.corners = corners
        gdf = gpd.GeoDataFrame({'index': ['bin']}, geometry=[Polygon(corners)],
                               index=pd.Index(['bin'], name='barcode'))
        adata = ad.AnnData(np.ones((1, 1)), obs=pd.DataFrame(index=['bin']),
                           var=pd.DataFrame(index=['gene']))
        self.sample = Bins(gdf, adata, pixels_per_bin=2, tile_scale_yfactor=2,
                           native_angle=89, bin_size=16)
        self.view = make_view(self.sample)
        self.zoom(20)

    def zoom(self, span):
        return show(self.view, 0, span, 0, span)

    def factor(self):
        factor = self.sample.get_aggregation_factor(self.view.rows.effective_resolution_level)
        self.assertEqual(self.view.bin_info()['factor'], factor)        # what the overlay shows
        return factor

    def use_clusters(self, *columns):
        """Bind cluster columns to the layer (colour and an extra width parameter)."""
        for column in columns:
            self.sample.adata.obs[column] = ['0']
        self.view.colour_by_categories(columns[0], {'0': '#ff0000'}, 'clusters')
        if len(columns) > 1:
            self.sample.adata.obs[columns[1]] = pd.Categorical(['0'])
            self.view.set_parameter('line_color', mode='data', column=columns[1])

    def test_native_dimensions_and_rotation_use_polygon_edges(self):
        info = self.view.header([])['bins']
        edge = self.corners[1] - self.corners[0]
        self.assertEqual(info['factor'], 1)
        np.testing.assert_allclose(info['width'], np.linalg.norm(edge))
        np.testing.assert_allclose(info['height'], np.linalg.norm(self.corners[2] - self.corners[1]))
        self.assertAlmostEqual(info['rotation'], math.degrees(math.atan2(edge[1], edge[0])))
        self.assertEqual(info['bin_size_um'], 16)
        self.assertIsNone(info['fixed'])
        self.assertIsNone(info['minimum'])
        self.assertEqual(info['factors'], [1, 2, 4, 8, 12])

    def test_zoom_updates_level_and_glyph_size(self):
        self.assertEqual(self.factor(), 1)
        self.assertEqual(self.view.header([])['level'], None)
        self.assertEqual(self.view.describe()['aggregation_factors'], [1, 2, 4, 8, 12])
        ops = self.zoom(10000)
        self.assertEqual(self.factor(), 12)
        header = self.view.header(ops)
        self.assertEqual(header['level'], 'agg:12')
        self.assertEqual(header['constants'],
                         dict(width=24.0, height=48.0, angle=0.0))   # rectangles at aggregated levels
        info = header['bins']
        self.assertEqual((info['factor'], info['width'], info['height'], info['rotation']), (12, 24, 48, 0))
        self.zoom(20)
        self.assertEqual(self.factor(), 1)
        self.assertEqual(loaded(self.view), ['bin'])

    def test_cluster_minimum_allows_coarser_zoom(self):
        self.use_clusters('leiden_4')
        self.zoom(10000)
        self.assertEqual(self.factor(), 12)
        self.assertEqual(self.view.rows.minimum_aggregation_factor, 4)
        self.assertEqual(self.view.bin_info()['minimum'], 4)
        self.assertIsNone(self.view.describe()['fixed_aggregation'])
        self.assertIsNone(self.view.bin_info()['fixed'])
        self.zoom(20)
        self.assertEqual(self.factor(), 4)                     # not below the analysis level
        self.view.colour_by_column('gene')
        self.assertIsNone(self.view.rows.minimum_aggregation_factor)
        self.zoom(20)
        self.assertEqual(self.factor(), 1)

    def test_fixed_level_survives_zoom_then_dynamic_uses_current_view(self):
        self.view.set_aggregation(4)
        self.zoom(20)
        self.assertEqual(self.view.rows.fixed_aggregation_factor, 4)
        self.assertEqual(self.view.describe()['fixed_aggregation'], 4)
        self.assertEqual(self.factor(), 4)
        self.assertTrue(all(str(name).endswith(':4') for name in loaded(self.view)))
        self.zoom(10000)
        self.assertEqual(self.factor(), 4)
        self.view.set_aggregation(None)
        self.zoom(10000)
        self.assertIsNone(self.view.rows.fixed_aggregation_factor)
        self.assertEqual(self.factor(), 12)

    def test_raw_mode_overrides_zoom_and_reloads_native_bins(self):
        self.zoom(10000)
        self.view.set_aggregation(4)
        self.view.set_aggregation(1)
        self.zoom(10000)
        self.assertEqual(self.view.rows.fixed_aggregation_factor, 1)
        self.assertIsNone(self.view.rows.effective_resolution_level)
        self.assertEqual(list(loaded(self.view)), ['bin'])
        self.assertEqual(self.view.rows.glyph_constants['width'], self.sample.native_width)
        np.testing.assert_allclose(self.view.bin_info()['width'],
                                   np.linalg.norm(self.corners[1] - self.corners[0]))

    def test_cluster_minimum_keeps_the_saved_fixed_level(self):
        self.view.set_aggregation(2)
        self.use_clusters('leiden_4')
        self.zoom(20)
        self.assertEqual(self.factor(), 4)
        self.assertEqual(self.view.rows.fixed_aggregation_factor, 2)
        with self.assertRaisesRegex(ValueError, 'at least 4'):        # a finer level is rejected
            self.view.set_aggregation(1)
        self.assertEqual(self.view.rows.fixed_aggregation_factor, 2)
        self.view.colour_by_column('gene')
        self.zoom(20)
        self.assertEqual(self.factor(), 2)

    def test_coarser_fixed_levels_and_dynamic_mode_remain_available(self):
        self.use_clusters('kmeans_2')
        for factor in (2, 4, 8, 12):
            self.view.set_aggregation(factor)
            self.zoom(20)
            self.assertEqual(self.factor(), factor)
        self.view.set_aggregation(None)
        self.zoom(20)
        self.assertIsNone(self.view.rows.fixed_aggregation_factor)
        self.assertEqual(self.factor(), 2)
        self.zoom(10000)
        self.assertEqual(self.factor(), 12)
        with self.assertRaisesRegex(ValueError, 'at least 2'):
            self.view.rows.set_aggregation_factor(1)

    def test_strictest_cluster_column_sets_minimum(self):
        for columns in (['leiden_2', 'kmeans_8'], ['kmeans_8', 'leiden_2']):
            view = self.view = make_view(self.sample)
            self.use_clusters(*columns)
            self.assertEqual(view.rows.minimum_aggregation_factor, 8)
            view.rows.set_fetch_context(bbox=[0, 20, 0, 20], resolution_level='agg:4')
            self.assertEqual(self.factor(), 8)

    def test_supported_fixed_levels_and_invalid_factor(self):
        for factor in self.sample.AGGREGATION_FACTORS:
            self.view.set_aggregation(factor)
            self.zoom(20)
            self.assertEqual(self.factor(), factor)
        with self.assertRaises(ValueError):
            self.view.set_aggregation(3)
        self.assertEqual(self.view.rows.fixed_aggregation_factor, 12)

    def test_empty_geometry(self):
        self.sample.gdf = self.sample.gdf.iloc[:0]
        info = make_view(self.sample).bin_info()
        self.assertIsNone(info['width'])
        self.assertIsNone(info['height'])
        self.assertIsNone(info['rotation'])

    def test_other_layers_have_no_bin_info(self):
        self.assertIsNone(make_view(kind='cells').header([])['bins'])


if __name__ == '__main__':
    unittest.main()
