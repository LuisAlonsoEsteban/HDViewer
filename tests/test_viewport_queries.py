"""Indexed queries retain the original viewport and aggregation semantics."""
import unittest

import hd_viewer.startup
import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Polygon, box

from hd_viewer.startup.data.image.transform import SpatialTransform
from hd_viewer.startup.sample import Bins, Cells
from hd_viewer.startup.sample._sample import Sample
from hd_viewer.startup.sample.utils.viewport import BlockViewportIndex


def bins_fixture(angle=0., reflected=False, side=12, shuffled=False, holes=False):
    row, col = np.indices((side, side))
    row, col = row.ravel(), col.ravel()
    if holes:
        keep = (row + col) % 5 != 0
        row, col = row[keep], col[keep]
    if shuffled:
        permutation = np.random.default_rng(42).permutation(len(row))
        row, col = row[permutation], col[permutation]
    theta = np.deg2rad(angle)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    basis = rotation @ np.array([[2., .4], [0., -3. if reflected else 3.]])
    matrix = np.eye(3)
    matrix[:2, :2] = basis
    matrix[:2, 2] = [30., 45.]
    identity = tuple(np.eye(3).ravel())
    transform = SpatialTransform(identity, identity, tuple(matrix.ravel()), 200, 1., 2.5, 'fixture')
    xy = transform.array_points(np.column_stack((col, row)))
    names = [f'b{i}' for i in range(len(row))]
    gdf = gpd.GeoDataFrame(dict(index=names, array_col=col, array_row=row),
                           index=names, geometry=transform.rectangles(xy))
    values = np.arange(len(row), dtype=float)
    values[::11] = np.nan
    adata = ad.AnnData(values[:, None], obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['gene']))
    adata.obs['fraction'] = np.arange(len(row), dtype=float) / max(1, len(row))
    adata.obs['kind'] = pd.Categorical(np.where(col % 2, 'a', 'b'))
    adata.uns['mean_obs_columns'] = ['fraction']
    return Bins(gdf, adata, bin_size=2, pixels_per_bin=2, spatial_transform=transform)


def old_grid_crop(sample, factor, bbox):
    grid = sample.aggregation_grid(factor)
    if bbox is None:
        return sample.gdf, grid.copy()
    left, right, bottom, top = bbox
    angle = sample.spatial_transform.angle if sample.spatial_transform else 0.
    width, height = sample.native_width * factor, sample.native_height * factor
    hx = (abs(np.cos(angle))*width + abs(np.sin(angle))*height)/2
    hy = (abs(np.sin(angle))*width + abs(np.cos(angle))*height)/2
    mask = ((grid.x + hx >= left) & (grid.x - hx <= right)
            & (grid.y + hy >= bottom) & (grid.y - hy <= top))
    return sample.gdf.loc[mask], grid.loc[mask].copy()


class NativeViewportTests(unittest.TestCase):
    def test_exact_predicate_order_and_boundary_contacts(self):
        # A bounding-box-only match must not include the hole in this polygon.
        ring = Polygon([(0, 0), (6, 0), (6, 6), (0, 6)],
                       holes=[[(1, 1), (5, 1), (5, 5), (1, 5)]])
        gdf = gpd.GeoDataFrame({'index': ['z', 'a', 'q', 'empty']}, index=['z', 'a', 'q', 'empty'],
                              geometry=[ring, box(10, 10, 12, 12), box(-3, -3, -1, -1), None])
        sample = Sample(gdf)
        for left, right, top, bottom in ((2, 3, 2, 3), (6, 7, 1, 2), (-10, 20, -10, 20),
                                        (100, 101, 100, 101), (7, -4, 8, -4),
                                        (None, 4, None, 4), (None, None, None, None), (6, 6, 1, 2)):
            with self.subTest(bounds=(left, right, top, bottom)):
                expected = gdf.cx[left:right, bottom:top].copy()
                actual = sample.crop_geom(left, right, top, bottom)
                pd.testing.assert_frame_equal(actual, expected)

    def test_native_bins_preserve_existing_padding_and_invalid_ranges(self):
        sample = bins_fixture(angle=27., reflected=True, shuffled=True)
        for bbox in ((25, 43, 20, 50), (30, 30, 45, 60), (30, 35, 60, 45), (None, 50, 20, 80)):
            left, right, top, bottom = bbox
            half = sample.bin_size / 2
            left, right = half if left is None else left, half if right is None else right
            top, bottom = half if top is None else top, half if bottom is None else bottom
            expected = (sample.gdf.iloc[:0] if left >= right or top >= bottom else
                        sample.gdf.cx[left-half:right+half, top-half:bottom+half].copy())
            pd.testing.assert_frame_equal(sample.crop_geom(*bbox), expected)

    def test_duplicate_labels_and_empty_layers(self):
        gdf = gpd.GeoDataFrame({'index': ['a', 'a']}, index=['a', 'a'],
                              geometry=[box(0, 0, 1, 1), box(3, 3, 4, 4)])
        sample = Sample(gdf)
        pd.testing.assert_frame_equal(sample.crop_geom(-1, 5, -1, 5), gdf)
        sample.gdf = gdf.iloc[:0].copy()
        pd.testing.assert_frame_equal(sample.crop_geom(-1, 5, -1, 5), sample.gdf)

    def test_spatial_index_tracks_geometry_edits(self):
        gdf = gpd.GeoDataFrame({'index': ['a']}, index=['a'], geometry=[box(0, 0, 1, 1)])
        sample = Sample(gdf)
        self.assertEqual(len(sample.crop_geom(-1, 2, -1, 2)), 1)
        sample.gdf.at['a', 'geometry'] = box(20, 20, 21, 21)
        self.assertTrue(sample.crop_geom(-1, 2, -1, 2).empty)
        self.assertEqual(len(sample.crop_geom(19, 22, 19, 22)), 1)

    def test_warm_query_reuses_tree_and_does_not_scan_all_intersections(self):
        sample = bins_fixture()
        sample.crop_geom(25, 40, 30, 50)
        tree = sample.gdf.sindex
        from unittest.mock import patch
        with patch.object(gpd.GeoDataFrame, 'intersects', side_effect=AssertionError('full scan')):
            sample.crop_geom(26, 41, 31, 51)
        self.assertIs(sample.gdf.sindex, tree)


class AggregatedViewportTests(unittest.TestCase):
    def assert_matches_scan(self, sample, factor, bbox):
        expected_gdf, expected_grid = old_grid_crop(sample, factor, bbox)
        actual_gdf, actual_grid = sample._crop_aggregation_grid(factor, bbox)
        pd.testing.assert_frame_equal(actual_gdf, expected_gdf)
        pd.testing.assert_frame_equal(actual_grid, expected_grid)

    def test_rotations_reflections_holes_and_native_order(self):
        rng = np.random.default_rng(7)
        for angle in (0., 27., 89.9, 179.):
            for reflected in (False, True):
                sample = bins_fixture(angle, reflected, shuffled=True, holes=True)
                for factor in (2, 4, 8):
                    for _ in range(12):
                        x, y = rng.uniform(-20, 100, 2)
                        self.assert_matches_scan(sample, factor, [x, x + 15, y, y + 20])
                    self.assert_matches_scan(sample, factor, [-1000, 1000, -1000, 1000])
                    self.assert_matches_scan(sample, factor, None)

    def test_exact_and_adjacent_float_boundaries(self):
        sample = bins_fixture(angle=27.)
        factor = 4
        grid = sample.aggregation_grid(factor)
        hx = (abs(np.cos(sample.spatial_transform.angle))*sample.native_width*factor
              + abs(np.sin(sample.spatial_transform.angle))*sample.native_height*factor)/2
        hy = (abs(np.sin(sample.spatial_transform.angle))*sample.native_width*factor
              + abs(np.cos(sample.spatial_transform.angle))*sample.native_height*factor)/2
        x, y = grid.iloc[0][['x', 'y']]
        for xedge in (x-hx, x+hx):
            for edge in (np.nextafter(xedge, -np.inf), xedge, np.nextafter(xedge, np.inf)):
                self.assert_matches_scan(sample, factor, [edge, edge+1, y-hy, y+hy])

    def test_complete_block_values_and_categorical_ties_match_original_fetch(self):
        # Run the original scan in the same aggregation pipeline as an oracle.
        from unittest.mock import patch
        for angle in (0., 27.):
            sample = bins_fixture(angle=angle, shuffled=True, holes=True)
            for factor in (2, 4):
                bbox = [29.9, 42.1, 45.1, 59.7]
                names = ['gene', 'fraction', 'kind']
                actual = sample.fetch_data(bbox=bbox, colnames=names, resolution_level=f'agg:{factor}')
                sample._multiscale_cache.clear()
                with patch.object(sample, '_crop_aggregation_grid',
                                  side_effect=lambda f, b: old_grid_crop(sample, f, b)):
                    expected = sample.fetch_data(bbox=bbox, colnames=names, resolution_level=f'agg:{factor}')
                pd.testing.assert_frame_equal(actual, expected)

    def test_index_has_one_point_per_block_and_is_reused(self):
        sample = bins_fixture()
        self.assert_matches_scan(sample, 4, [30, 40, 45, 55])
        index = sample._viewport_indexes[4]
        self.assertEqual(len(index.tree), len(sample.aggregation_grid(4).drop_duplicates(['grid_x', 'grid_y'])))
        self.assertLess(len(index.tree), len(sample.gdf))
        self.assert_matches_scan(sample, 4, [31, 41, 46, 56])
        self.assertIs(sample._viewport_indexes[4], index)

    def test_replacing_geometry_frame_or_grid_rebuilds_index(self):
        sample = bins_fixture()
        self.assert_matches_scan(sample, 4, [30, 40, 45, 55])
        old = sample._viewport_indexes[4]
        sample.gdf = sample.gdf.iloc[::-1].copy()
        self.assert_matches_scan(sample, 4, [30, 40, 45, 55])
        self.assertIsNot(sample._viewport_indexes[4], old)
        old = sample._viewport_indexes[4]
        sample._grid_cache[4] = sample._grid_cache[4].copy()
        self.assert_matches_scan(sample, 4, [30, 40, 45, 55])
        self.assertIsNot(sample._viewport_indexes[4], old)

    def test_empty_grid(self):
        grid = pd.DataFrame({name: np.array([], dtype=float) for name in ('x', 'y', 'grid_x', 'grid_y')})
        index = BlockViewportIndex(grid)
        self.assertEqual(index.query(0, 1, 0, 1).tolist(), [])

    def test_in_memory_bins_without_transform(self):
        sample = bins_fixture(angle=27., shuffled=True, holes=True)
        sample.spatial_transform = None
        for factor in (2, 4, 8):
            for bbox in ([30, 45, 45, 60], [-1000, 1000, -1000, 1000], [300, 301, 300, 301]):
                self.assert_matches_scan(sample, factor, bbox)

    def test_empty_fetch_and_unsequenced_members_match_scan(self):
        from unittest.mock import patch
        sample = bins_fixture(shuffled=True)
        sample.adata = sample.adata[::3].copy()
        for bbox in ([30, 42, 45, 60], [300, 301, 300, 301]):
            actual = sample.fetch_data(bbox=bbox, colnames=['gene', 'kind'], resolution_level='agg:4')
            sample._multiscale_cache.clear()
            with patch.object(sample, '_crop_aggregation_grid',
                              side_effect=lambda f, b: old_grid_crop(sample, f, b)):
                expected = sample.fetch_data(bbox=bbox, colnames=['gene', 'kind'], resolution_level='agg:4')
            pd.testing.assert_frame_equal(actual, expected)


if __name__ == '__main__':
    unittest.main()
