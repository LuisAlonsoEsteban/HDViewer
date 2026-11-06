"""Native barcode search, calibration, and radius ROI lifecycle."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import shapely

from hd_viewer import settings
from tests.test_roi_manager import Api, make_layer


def barcode(identifier, radius=.2, **extra):
    return dict(kind='barcode', barcode=identifier, radius=radius, **extra)


def bounds(preview):
    return preview['bounds']


class BarcodeROITests(unittest.TestCase):
    def test_prefix_suggestions_are_bounded_and_server_index_is_reused(self):
        with patch.object(settings, 'ROI_BARCODE_MAX_COMPLETIONS', 3):
            api = Api(self)
            self.assertIsNone(api.view.rois._barcode_cache)            # built only when first used
            self.assertEqual(api.ok(api.get('barcodes', params=dict(q='b'))), ['b0_0', 'b0_1', 'b0_2'])
            cached = api.view.rois._barcode_cache
            self.assertEqual(api.ok(api.get('barcodes', params=dict(q='b3_'))), ['b3_0', 'b3_1', 'b3_2'])
            self.assertIs(api.view.rois._barcode_cache, cached)
            self.assertEqual(api.ok(api.get('barcodes', params=dict(q=''))), [])
            self.assertEqual(api.ok(api.get('barcodes', params=dict(q='  '))), [])
        # Other layers have no native bin index to search; map seeds remain available
        spots = Api(self, 'spots')
        self.assertEqual(spots.ok(spots.get('barcodes', params=dict(q='s'))), [])
        self.assertEqual(api.roi(kind='point', centre=[.5, .5], radius=.2)['count'], 1)

    def test_search_uses_native_bin_outside_loaded_view_and_aggregate_display(self):
        api = Api(self)
        rows = api.view.rows
        rows.set_fetch_context(bbox=[0, 1, 0, 1], resolution_level='agg:4')
        rows.fetch_data()
        self.assertNotIn('b3_3', list(rows.loaded_obs_names))
        preview = api.ok(api.post('rois/preview', json=barcode('b3_3')))
        np.testing.assert_allclose(bounds(preview), [3.3, 3.3, 3.7, 3.7])
        self.assertIsNone(preview['seeds'])
        geometry = api.view.rois.seed_radius_geometry(3.5, 3.5, .2)
        self.assertEqual(int(api.view.rois.compute_membership(geometry).sum()), 1)
        self.assertEqual(api.roi(**barcode('b3_3'))['count'], 1)

    def test_saved_roi_records_barcode_and_radius_changes_do_not_edit_saved_region(self):
        api = Api(self)
        roi = api.roi(**barcode('b3_3', name='Barcode region'))
        saved = api.view.rois.rois[0]
        geometry = saved.geometry
        self.assertEqual((roi['name'], roi['seed_barcode'], roi['count']), ('Barcode region', 'b3_3', 1))
        self.assertEqual(saved.seed_barcode, 'b3_3')
        larger = api.ok(api.post('rois/preview', json=barcode('b3_3', .8)))
        self.assertGreater(shapely.box(*bounds(larger)).area, shapely.box(*geometry.bounds).area)
        self.assertIs(saved.geometry, geometry)
        self.assertEqual(api.ok(api.get('rois'))['rois'][0]['seed_barcode'], 'b3_3')

    def test_invalid_or_partial_barcode_cannot_be_previewed_or_saved(self):
        api = Api(self)
        for name, message in (('b3_', 'not found'), ('', 'Enter a native barcode')):
            with self.subTest(barcode=name):
                for path in ('rois/preview', 'rois'):
                    response = api.post(path, json=barcode(name))
                    self.assertEqual(response.status_code, 422)
                    self.assertIn(message, response.json()['detail'])
        self.assertEqual(api.ok(api.get('rois'))['rois'], [])

    def test_changing_radius_revalidates(self):
        api = Api(self)
        response = api.post('rois/preview', json=barcode('b2_2', 0))
        self.assertEqual(response.status_code, 422)
        self.assertIn('positive', response.json()['detail'])
        self.assertEqual(api.post('rois/preview', json=dict(kind='barcode', barcode='b2_2')).status_code, 422)
        np.testing.assert_allclose(bounds(api.ok(api.post('rois/preview', json=barcode('b2_2', .5)))),
                                   [2., 2., 3., 3.])

    def test_physical_calibration_matches_manual_seed_radius(self):
        api = Api(self)
        rois = api.view.rois
        rois.coordinate_scale = SimpleNamespace(data=dict(px_x=[2.], px_y=[.5], um_per_px=[.5]))
        self.assertEqual(api.ok(api.get('rois'))['unit'], 'µm')
        preview = api.ok(api.post('rois/preview', json=barcode('b1_1', 1.)))
        expected = rois.seed_radius_geometry(1.5, 1.5, 1., *rois.radius_scale()[0])
        np.testing.assert_allclose(bounds(preview), expected.bounds)
        np.testing.assert_allclose(bounds(preview), [.5, -2.5, 2.5, 5.5])
        api.roi(**barcode('b1_1', 1.))
        self.assertTrue(rois.rois[0].geometry.equals_exact(expected, 1e-12))
        # A map click at the same centre gives the same region
        api.roi(kind='point', centre=[1.5, 1.5], radius=1.)
        self.assertTrue(rois.rois[1].geometry.equals_exact(expected, 1e-12))

    def test_transformed_image_uses_the_same_seed_and_calibration_as_a_map_click(self):
        from hd_viewer.startup.data.image.transform import SpatialTransform
        api = Api(self)
        sample, rois = api.view.sample, api.view.rois
        matrix = np.array([[.5, 0., 10.], [0., -.5, 20.], [0., 0., 1.]])
        transform = SpatialTransform(tuple(matrix.ravel()), tuple(np.linalg.inv(matrix).ravel()),
                                     tuple(matrix.ravel()), 40, .5, .5, 'test')
        sample.gdf = sample.gdf.copy()
        sample.gdf.geometry = shapely.transform(sample.gdf.geometry.to_numpy(), transform.points)
        sample.spatial_transform = transform
        rois.coordinate_scale = SimpleNamespace(data=dict(px_x=[2.], px_y=[2.], um_per_px=[.5]))
        np.testing.assert_allclose(rois.barcode_seed_center('b1_1'), [10.75, 19.25])
        preview = api.ok(api.post('rois/preview', json=barcode('b1_1', 1.)))
        np.testing.assert_allclose(bounds(preview), [9.75, 18.25, 11.75, 20.25])
        click = api.ok(api.post('rois/preview', json=dict(kind='point', centre=[10.75, 19.25], radius=1.)))
        np.testing.assert_allclose(bounds(click), bounds(preview))

    def test_exact_identifier_validation_and_cache_refresh(self):
        view = make_layer()
        rois, sample = view.rois, view.sample
        self.assertEqual(rois.barcode_seed_center(' b0_0 '), (.5, .5))
        cached = rois._barcode_cache
        sample.gdf = sample.gdf.iloc[::-1].copy()
        self.assertEqual(rois.barcode_seed_center('b0_0'), (.5, .5))
        self.assertIsNot(rois._barcode_cache, cached)
        sample.gdf['index'] = ['replacement'] * len(sample.gdf)
        with self.assertRaisesRegex(ValueError, 'more than one'):
            rois.barcode_seed_center('replacement')
        with self.assertRaisesRegex(ValueError, 'not found'):
            rois.barcode_seed_center('b0_0')

    def test_native_numeric_identifiers_and_unsequenced_bins_can_supply_seeds(self):
        api = Api(self)
        rois, sample = api.view.rois, api.view.sample
        np.testing.assert_allclose(rois.barcode_seed_center('not_sequenced0'), [1.3, 1.3])
        self.assertEqual(api.roi(**barcode('not_sequenced0', .05))['count'], 1)   # only b1_1, around it
        sample.gdf['index'] = np.arange(len(sample.gdf))
        self.assertIn('10', api.ok(api.get('barcodes', params=dict(q='1'))))
        self.assertEqual(rois.barcode_seed_center('0'), (.5, .5))

    def test_missing_geometry_is_handled_explicitly(self):
        api = Api(self)
        api.view.sample.gdf.at['b0_0', 'geometry'] = None
        response = api.post('rois/preview', json=barcode('b0_0'))
        self.assertEqual(response.status_code, 422)
        self.assertIn('valid seed centre', response.json()['detail'])


if __name__ == '__main__':
    unittest.main()
