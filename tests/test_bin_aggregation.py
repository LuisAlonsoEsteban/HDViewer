import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from bokeh_server.manager import DynamicSourceManager
from src.sample import Bins


class BinAggregationGeometryTests(unittest.TestCase):
    def make_sample(self, scale_x=1.0, scale_y=2.0):
        width, height = 2 * scale_x, 2 * scale_y
        ids = [str(i) for i in range(64)]
        polygons = [box(x * width, y * height, (x + 1) * width, (y + 1) * height)
                    for y in range(8) for x in range(8)]
        gdf = gpd.GeoDataFrame({'index': ids}, geometry=polygons,
                               index=pd.Index(ids, name='barcode'))
        adata = ad.AnnData(np.ones((64, 1)), obs=pd.DataFrame(index=ids),
                          var=pd.DataFrame(index=['gene']))
        return Bins(gdf, adata, pixels_per_bin=2, native_angle=1.55566,
                    tile_scale_xfactor=scale_x, tile_scale_yfactor=scale_y)

    def test_aggregated_rectangles_tile_both_axes(self):
        for scale_x, scale_y in ((1, 2), (2, 1), (1, 1)):
            sample = self.make_sample(scale_x, scale_y)
            for factor in (2, 4):
                data = sample.fetch_data(resolution_level=f'agg:{factor}', colnames=['gene'])
                self.assertTrue((data['angle'] == 0).all())
                self.assertTrue((data['width'] == sample.native_width * factor).all())
                self.assertTrue((data['height'] == sample.native_height * factor).all())
                np.testing.assert_allclose(np.diff(sorted(data.x.unique())), data.width.iloc[0])
                np.testing.assert_allclose(np.diff(sorted(data.y.unique())), data.height.iloc[0])
                self.assertEqual(data.gene.sum(), 64)
                cached = sample.fetch_data(resolution_level=f'agg:{factor}', colnames=['gene'])
                pd.testing.assert_frame_equal(data, cached)

    def test_native_rotation_and_resolution_switching(self):
        sample = self.make_sample()
        manager = DynamicSourceManager(sample)
        for level, expected_angle in ((None, sample.native_angle), ('agg:2', 0),
                                      ('agg:4', 0), (None, sample.native_angle)):
            manager.set_fetch_context(resolution_level=level)
            manager.clear_sources()
            manager.add_data()
            self.assertTrue(manager.source.data['angle'])
            self.assertTrue(all(a == expected_angle for a in manager.source.data['angle']))

    def test_renderer_uses_native_and_aggregated_corners(self):
        from bokeh.document import Document
        from bokeh.plotting import figure
        from bokeh.models import Patches
        from bokeh_server.components.glyphs import Squares

        sample = self.make_sample()
        plot = figure()
        component = Squares(plot, sample, lazy_loading=True, alias='Bin', hover_color='green')
        component.manager.clear_sources()
        component.add_glyph()
        component.manager.add_data()
        glyph = component.glyph_renderer.glyph
        self.assertIsInstance(glyph, Patches)
        self.assertEqual(glyph.xs, 'xs')
        self.assertEqual(glyph.ys, 'ys')
        component.manager.set_fetch_context(resolution_level='agg:2')
        component.manager.add_data()
        doc = Document()
        doc.add_root(plot)
        doc.to_json()

    def test_empty_crop_keeps_angle_schema(self):
        sample = self.make_sample()
        data = sample.fetch_data(bbox=[100, 110, 100, 110], resolution_level='agg:2')
        self.assertTrue(data.empty)
        self.assertIn('angle', data.columns)

class NativeBinFootprintTests(unittest.TestCase):
    def test_rotated_reflected_grids_with_unequal_axis_scaling(self):
        from src.utils.bins2spots import Bins2Spots
        from shapely.geometry import Polygon

        for angle in (0.0, np.deg2rad(179), np.deg2rad(25)):
            for scale_x, scale_y in ((1, 2), (2, 1), (1, 1)):
                with self.subTest(angle=angle, scales=(scale_x, scale_y)):
                    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                                         [np.sin(angle), np.cos(angle)]])
                    # Include a reflection, as in the example tissue positions.
                    basis = np.diag([2.0, -2.0]) @ rotation.T
                    rows, cols = np.indices((4, 4))
                    centers = np.column_stack((cols.ravel(), rows.ravel())) @ basis + [30, 40]
                    ids = [str(i) for i in range(16)]
                    positions = pd.DataFrame({
                        'array_row': rows.ravel(), 'array_col': cols.ravel(),
                        'pxl_col_in_fullres': centers[:, 0],
                        'pxl_row_in_fullres': centers[:, 1], 'index': ids,
                    }, index=pd.Index(ids, name='barcode'))
                    b2s = Bins2Spots(None, positions, bins_size=2, pixels_per_bin=2,
                                    tile_scale_xfactor=scale_x, tile_scale_yfactor=scale_y)
                    polygons = b2s.bin_squares
                    expected_centers = centers * [scale_x, scale_y]
                    np.testing.assert_allclose(polygons.geometry.centroid.x, expected_centers[:, 0])
                    np.testing.assert_allclose(polygons.geometry.centroid.y, expected_centers[:, 1])
                    np.testing.assert_allclose(polygons.geometry.area, 4 * scale_x * scale_y)
                    first = polygons.geometry.iloc[0]
                    for neighbor in (1, 4):
                        other = polygons.geometry.iloc[neighbor]
                        self.assertLess(first.distance(other), 1e-10)
                        self.assertLess(first.intersection(other).area, 1e-10)
                    adata = ad.AnnData(np.ones((16, 1)), obs=pd.DataFrame(index=ids),
                                      var=pd.DataFrame(index=['gene']))
                    sample = Bins(polygons, adata, pixels_per_bin=2,
                                  tile_scale_xfactor=scale_x, tile_scale_yfactor=scale_y)
                    native = sample.fetch_data(colnames=['gene'])
                    for index, item in native.iterrows():
                        rendered = Polygon(zip(item['xs'], item['ys']))
                        self.assertTrue(rendered.equals(polygons.loc[index].geometry))
                    manager = DynamicSourceManager(sample)
                    for level in (None, 'agg:2', None):
                        manager.set_fetch_context(resolution_level=level)
                        manager.clear_sources()
                        manager.add_data()
                        self.assertEqual(len(manager.source.data['xs']), len(manager.source.data['x']))
                        self.assertTrue(all(len(xs) == 5 for xs in manager.source.data['xs']))


if __name__ == '__main__':
    unittest.main()
