"""Native barcode search, calibration, and radius ROI lifecycle."""
import unittest
from unittest.mock import patch

import hd_viewer.startup  # Match the app's import order.
from hd_viewer import settings
import numpy as np
import shapely
from bokeh.events import ButtonClick
from bokeh.models import ColumnDataSource

from tests.test_roi_manager import make_layer, wait_for


def barcode_mode(layer):
    manager = layer.roi_manager
    manager.method.value = 'radius'
    manager.radius.value = .2
    manager.seed_source.value = 'barcode'
    return manager


class BarcodeROITests(unittest.TestCase):
    def test_control_visibility_and_map_mode_remain_available(self):
        layer = make_layer()
        manager = layer.roi_manager
        self.assertEqual(manager.seed_barcode.completions, [])
        self.assertIsNone(manager._barcode_cache)
        self.assertFalse(manager.seed_barcode.visible)
        barcode_mode(layer)
        self.assertTrue(manager.seed_barcode.visible)
        self.assertTrue(manager.seed_hint.visible)
        self.assertFalse(manager.seed_group.visible)
        self.assertFalse(manager.group_region.visible)
        self.assertEqual(manager.draw_button.label, 'Preview barcode')
        manager.seed_source.value = 'map'
        self.assertFalse(manager.seed_barcode.visible)
        self.assertEqual(manager.draw_button.label, 'Draw')
        manager.begin()
        self.assertIs(manager.drawing.owner, manager)
        manager.cancel()

    def test_prefix_suggestions_are_bounded_and_server_index_is_reused(self):
        with patch.object(settings, 'ROI_BARCODE_MAX_COMPLETIONS', 3):
            layer = make_layer()
            manager = barcode_mode(layer)
            manager.seed_barcode.value_input = 'b'
            self.assertEqual(manager.seed_barcode.completions, ['b0_0', 'b0_1', 'b0_2'])
            self.assertEqual(manager.seed_barcode.max_completions, 3)
            cached = manager._barcode_cache
            manager.seed_barcode.value_input = 'b3_'
            self.assertEqual(manager.seed_barcode.completions, ['b3_0', 'b3_1', 'b3_2'])
            self.assertIs(manager._barcode_cache, cached)
            manager.seed_barcode.value_input = ''
            self.assertEqual(manager.seed_barcode.completions, [])

    def test_search_uses_native_bin_outside_loaded_view_and_aggregate_display(self):
        layer = make_layer()
        layer.manager.set_fetch_context(bbox=[0, 1, 0, 1], resolution_level='agg:4')
        layer.manager.fetch_data()
        self.assertNotIn('b3_3', list(layer.manager.loaded_obs_names))
        manager = barcode_mode(layer)
        manager.seed_barcode.value = 'b3_3'
        np.testing.assert_allclose(manager.draft.bounds, [3.3, 3.3, 3.7, 3.7])
        self.assertEqual(manager._draft_seed_barcode, 'b3_3')
        self.assertIsNone(manager.drawing.owner)
        self.assertFalse(manager.save_button.disabled)
        self.assertLess(layer.figure.x_range.start, 3.3)
        self.assertGreater(layer.figure.x_range.end, 3.7)
        self.assertEqual(int(manager.compute_membership(manager.draft).sum()), 1)

    def test_saved_roi_records_barcode_and_radius_changes_do_not_edit_saved_region(self):
        manager = barcode_mode(make_layer())
        manager.seed_barcode.value = 'b3_3'
        geometry = manager.draft
        manager.name_input.value = 'Barcode region'
        manager.save_roi()
        wait_for(lambda: not manager.logger.running)
        roi = manager.rois[0]
        self.assertEqual(roi.seed_barcode, 'b3_3')
        self.assertEqual(roi.count, 1)
        self.assertIs(roi.geometry, geometry)
        self.assertIn('b3_3', manager.roi_list.children[0].children[0].text)
        manager.radius.value = .8
        self.assertGreater(manager.draft.area, roi.geometry.area)
        self.assertIs(roi.geometry, geometry)
        self.assertEqual(roi.seed_barcode, 'b3_3')

    def test_invalid_or_partial_barcode_clears_draft_and_cannot_be_saved(self):
        manager = barcode_mode(make_layer())
        manager.seed_barcode.value = 'b0_0'
        manager.seed_barcode.value_input = 'b3_'
        self.assertIsNone(manager.draft)
        self.assertTrue(manager.save_button.disabled)
        manager.seed_barcode.value = 'b3_'
        self.assertIsNone(manager.draft)
        self.assertIn('not found', manager.logger.div.text)
        manager.save_roi()
        self.assertEqual(manager.rois, [])
        manager.seed_barcode.value = ''
        manager.begin()
        self.assertIn('Choose a barcode', manager.logger.div.text)

    def test_changing_radius_revalidates_and_preview_button_refocuses(self):
        layer = make_layer()
        manager = barcode_mode(layer)
        manager.seed_barcode.value = 'b2_2'
        manager.radius.value = 0
        self.assertIsNone(manager.draft)
        self.assertTrue(manager.save_button.disabled)
        manager.radius.value = .5
        np.testing.assert_allclose(manager.draft.bounds, [2., 2., 3., 3.])
        layer.figure.x_range.update(start=100, end=200)
        manager.draw_button._trigger_event(ButtonClick(manager.draw_button))
        self.assertLess(layer.figure.x_range.start, 2.)
        self.assertGreater(layer.figure.x_range.end, 3.)

    def test_physical_calibration_matches_manual_seed_radius(self):
        manager = barcode_mode(make_layer())
        manager.coordinate_scale = ColumnDataSource(dict(px_x=[2.], px_y=[.5], um_per_px=[.5]))
        manager.radius.value = 1.
        manager.seed_barcode.value = 'b1_1'
        expected = manager.seed_radius_geometry(1.5, 1.5, 1., *manager.radius_scale()[0])
        self.assertTrue(manager.draft.equals_exact(expected, 1e-12))
        np.testing.assert_allclose(manager.draft.bounds, [.5, -2.5, 2.5, 5.5])

    def test_transformed_image_uses_the_same_seed_and_calibration_as_a_map_click(self):
        from hd_viewer.startup.data.image.transform import SpatialTransform
        layer = make_layer()
        matrix = np.array([[.5, 0., 10.], [0., -.5, 20.], [0., 0., 1.]])
        transform = SpatialTransform(tuple(matrix.ravel()), tuple(np.linalg.inv(matrix).ravel()),
                                     tuple(matrix.ravel()), 40, .5, .5, 'test')
        layer.sample.gdf = layer.sample.gdf.copy()
        layer.sample.gdf.geometry = shapely.transform(layer.sample.gdf.geometry.to_numpy(), transform.points)
        layer.sample.spatial_transform = transform
        manager = barcode_mode(layer)
        manager.coordinate_scale = ColumnDataSource(dict(px_x=[2.], px_y=[2.], um_per_px=[.5]))
        manager.radius.value = 1.
        manager.seed_barcode.value = 'b1_1'
        np.testing.assert_allclose(manager.barcode_seed_center('b1_1'), [10.75, 19.25])
        np.testing.assert_allclose(manager.draft.bounds, [9.75, 18.25, 11.75, 20.25])

    def test_exact_identifier_validation_and_cache_refresh(self):
        layer = make_layer()
        manager = barcode_mode(layer)
        self.assertEqual(manager.barcode_seed_center(' b0_0 '), (.5, .5))
        cached = manager._barcode_cache
        layer.sample.gdf = layer.sample.gdf.iloc[::-1].copy()
        self.assertEqual(manager.barcode_seed_center('b0_0'), (.5, .5))
        self.assertIsNot(manager._barcode_cache, cached)
        layer.sample.gdf['index'] = ['replacement'] * len(layer.sample.gdf)
        with self.assertRaisesRegex(ValueError, 'more than one'):
            manager.barcode_seed_center('replacement')
        with self.assertRaisesRegex(ValueError, 'not found'):
            manager.barcode_seed_center('b0_0')

    def test_native_numeric_identifiers_and_unsequenced_bins_can_supply_seeds(self):
        layer = make_layer()
        manager = barcode_mode(layer)
        np.testing.assert_allclose(manager.barcode_seed_center('not_sequenced'), [1.3, 1.3])
        layer.sample.gdf['index'] = np.arange(len(layer.sample.gdf))
        manager.seed_barcode.value_input = '1'
        self.assertIn('10', manager.seed_barcode.completions)
        self.assertEqual(manager.barcode_seed_center('0'), (.5, .5))

    def test_missing_geometry_is_handled_explicitly(self):
        layer = make_layer()
        manager = barcode_mode(layer)
        layer.sample.gdf.at['b0_0', 'geometry'] = None
        manager.seed_barcode.value = 'b0_0'
        self.assertIsNone(manager.draft)
        self.assertIn('valid seed centre', manager.logger.div.text)


if __name__ == '__main__':
    unittest.main()
