import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from hd_viewer.core.layer_state import LayerRows
from hd_viewer.startup.sample import Bins


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
        rows = LayerRows(sample)
        for level, expected_angle in ((None, 0), ('agg:2', 0), ('agg:4', 0), (None, 0)):
            rows.set_fetch_context(resolution_level=level)
            rows.clear_sources()
            rows.add_data()
            # Constant per level: a glyph value, not a per-bin column
            self.assertNotIn('angle', rows.sink.data)
            self.assertEqual(rows.glyph_constants['angle'], expected_angle)
            factor = sample.get_aggregation_factor(level)
            self.assertEqual(rows.glyph_constants['width'], sample.native_width * factor)

    def test_rows_are_rectangle_centres_without_polygon_columns(self):
        """Bins are drawn as rectangles: the source holds centres only; width and
        height are per-level glyph constants that follow the aggregation level."""
        from tests.fixtures import make_view, show
        sample = self.make_sample()
        view = make_view(sample)
        show(view, -1, 17, -1, 33, width=800, height=800)
        data = view.rows.sink.data
        self.assertGreater(len(data['x']), 0)
        for name in ('width', 'height', 'angle', 'xs', 'ys'):
            self.assertNotIn(name, data)
        self.assertEqual(view.rows.glyph_constants['width'], sample.native_width)
        self.assertEqual(view.rows.glyph_constants['height'], sample.native_height)
        view.set_aggregation(2)
        show(view, -1, 17, -1, 33, width=800, height=800)
        self.assertEqual(view.rows.glyph_constants['width'], sample.native_width * 2)
        self.assertEqual(view.rows.glyph_constants['height'], sample.native_height * 2)
        header = view.header([])
        self.assertEqual(header['constants']['width'], sample.native_width * 2)
        self.assertEqual(header['level'], 'agg:2')

    def test_empty_crop_keeps_angle_schema(self):
        sample = self.make_sample()
        data = sample.fetch_data(bbox=[100, 110, 100, 110], resolution_level='agg:2')
        self.assertTrue(data.empty)
        self.assertIn('angle', data.columns)

class NativeBinFootprintTests(unittest.TestCase):
    def test_rotated_reflected_grids_with_unequal_axis_scaling(self):
        from hd_viewer.startup.operations.bins2spots import Bins2Spots
        from shapely.geometry import box
        from shapely.affinity import rotate
        from hd_viewer.startup.data.image.transform import SpatialTransform

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
                    transform = SpatialTransform.fit(positions, (100, 120, 3))
                    b2s = Bins2Spots(None, positions, bins_size=2, pixels_per_bin=2,
                                    spatial_transform=transform,
                                    tile_scale_xfactor=scale_x, tile_scale_yfactor=scale_y)
                    polygons = b2s.bin_squares
                    expected_centers = transform.points(centers)
                    np.testing.assert_allclose(polygons.geometry.centroid.x, expected_centers[:, 0])
                    np.testing.assert_allclose(polygons.geometry.centroid.y, expected_centers[:, 1])
                    np.testing.assert_allclose(polygons.geometry.area, transform.pitch ** 2)
                    first = polygons.geometry.iloc[0]
                    for neighbor in (1, 4):
                        other = polygons.geometry.iloc[neighbor]
                        self.assertLess(first.distance(other), 1e-10)
                        self.assertLess(first.intersection(other).area, 1e-10)
                    adata = ad.AnnData(np.ones((16, 1)), obs=pd.DataFrame(index=ids),
                                      var=pd.DataFrame(index=['gene']))
                    sample = Bins(polygons, adata, pixels_per_bin=2, spatial_transform=transform,
                                  tile_scale_xfactor=scale_x, tile_scale_yfactor=scale_y)
                    native = sample.fetch_data(colnames=['gene'])
                    for index, item in native.iterrows():
                        rendered = box(item.x-item.width/2, item.y-item.height/2,
                                       item.x+item.width/2, item.y+item.height/2)
                        rendered = rotate(rendered, item.angle, origin=(item.x,item.y), use_radians=True)
                        self.assertLess(rendered.symmetric_difference(polygons.loc[index].geometry).area, 1e-10)
                    rows = LayerRows(sample)
                    for level in (None, 'agg:2', None):
                        rows.set_fetch_context(resolution_level=level)
                        rows.clear_sources()
                        rows.add_data()
                        self.assertNotIn('xs', rows.sink.data)
                        self.assertNotIn('ys', rows.sink.data)


if __name__ == '__main__':
    unittest.main()
