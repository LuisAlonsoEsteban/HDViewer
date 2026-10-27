import threading
import unittest

from bokeh.document import Document
from bokeh.events import ButtonClick
from shapely.geometry import Polygon, box
from shapely import union_all
import numpy as np
import anndata as ad
import geopandas as gpd
import pandas as pd

from src.sample import Bins

from bokeh_server.components.utils.roi import ROI
from tests.test_roi_manager import ImmediateDoc, make_layer, wait_for


def create_seed_group(layer, names, name='Seeds'):
    return layer.groups_manager.add_mask_group(name, layer.sample.adata.obs_names.isin(names),
                                               layer.roi_manager.doc)


def fire_preview(manager):
    if manager._preview_timeout is not None:
        manager._preview_timeout.callback()
    wait_for(lambda: not manager.logger.running and manager._pending_preview is None)


class GroupRadiusROITests(unittest.TestCase):
    def test_large_group_uses_every_seed_across_union_batches(self):
        layer = make_layer()
        n = 23
        names = [f'b{x}_{y}' for y in range(n) for x in range(n)]
        gdf = gpd.GeoDataFrame({'index': names}, index=names,
                               geometry=[box(x, y, x + 1, y + 1) for y in range(n) for x in range(n)])
        data = ad.AnnData(np.ones((len(names), 1)), obs=pd.DataFrame(index=names),
                          var=pd.DataFrame(index=['g0']))
        data.obs['grp_All'] = True
        layer.manager.sample = Bins(gdf, data, pixels_per_bin=1, bin_size_um=2)
        shape, count = layer.roi_manager.compute_group_radius('All', .8, (1., 1.))
        self.assertEqual(count, 529)
        self.assertTrue(shape.is_valid)
        np.testing.assert_allclose(shape.bounds, [-.3, -.3, 23.3, 23.3])
        self.assertTrue(layer.roi_manager.compute_membership(shape).all())

    def test_all_group_bins_create_one_roi_even_outside_viewport_and_at_aggregate_level(self):
        layer = make_layer()
        manager = layer.roi_manager
        group = create_seed_group(layer, ['b0_0', 'b3_3'])
        # Hide the group and load an aggregated viewport away from one seed.
        layer.groups_manager._active_clusters.clear()
        layer.groups_manager._apply_visibility_mask()
        layer.manager.set_fetch_context(bbox=(0, 1, 0, 1), resolution_level='agg:2')
        layer.manager.fetch_data()
        manager.radius.value = .2
        manager.method.value = 'radius'
        manager.seed_source.value = 'group'
        manager.seed_group.value = group
        self.assertIsNotNone(manager._preview_timeout)  # Automatic, debounced preview.
        fire_preview(manager)
        self.assertEqual(manager._draft_seed_info, (group, 2))
        self.assertEqual(manager.draft.geom_type, 'MultiPolygon')
        self.assertEqual(len(manager.draft.geoms), 2)
        expected = union_all([manager.radius_geometry(.5, .5, .2), manager.radius_geometry(3.5, 3.5, .2)])
        self.assertLess(manager.draft.symmetric_difference(expected).area, 1e-10)
        self.assertEqual(len(manager.source.data['xs']), 1)  # One ROI, with two pieces.
        self.assertEqual(len(manager.source.data['xs'][0]), 2)
        self.assertEqual(manager.get_rois(), [])
        self.assertEqual(layer.get_selected_obs_names(), set())
        manager.name_input.value = 'Group neighbourhood'
        manager.save_button._trigger_event(ButtonClick(manager.save_button))
        wait_for(lambda: not manager.logger.running)
        self.assertEqual(len(manager.get_rois()), 1)
        roi = manager.get_rois()[0]
        self.assertEqual(roi.seed_count, 2)
        self.assertEqual(roi.seed_group, group)
        self.assertEqual(roi.count, 2)
        created = manager.create_group(roi)
        np.testing.assert_array_equal(layer.sample.adata.obs[f'grp_{created}'], roi.mask)
        manager.remove_roi(roi)
        self.assertIn(group, layer.groups_manager._clusters)
        self.assertNotIn(created, layer.groups_manager._clusters)

    def test_radius_changes_rebuild_preview_without_modifying_saved_rois(self):
        layer = make_layer()
        manager = layer.roi_manager
        group = create_seed_group(layer, ['b0_0', 'b1_0'])
        manager.method.value = 'radius'
        manager.seed_source.value = 'group'
        manager.radius.value = .2
        manager.seed_group.value = group
        fire_preview(manager)
        self.assertEqual(manager.draft.geom_type, 'MultiPolygon')
        manager.save_roi()
        wait_for(lambda: not manager.logger.running)
        roi = manager.get_rois()[0]
        saved_area = roi.geometry.area
        manager.radius.value = 1.
        fire_preview(manager)
        self.assertEqual(manager.draft.geom_type, 'Polygon')
        self.assertGreater(manager.draft.area, saved_area)
        self.assertEqual(roi.geometry.area, saved_area)
        self.assertEqual(len(manager.get_rois()), 1)
        # Membership changes of the selected group update the draft too.
        layer.sample.adata.obs[f'grp_{group}'] = layer.sample.adata.obs_names == 'b3_3'
        layer.groups_manager._notify_group_change()
        fire_preview(manager)
        self.assertEqual(manager._draft_seed_info, (group, 1))
        np.testing.assert_allclose(manager.draft.centroid.coords[0], [3.5, 3.5])

    def test_physical_scale_and_holes_are_preserved(self):
        layer = make_layer()
        manager = layer.roi_manager
        group = create_seed_group(layer, ['b0_0'])
        geometry, count = manager.compute_group_radius(group, 2., (2., .5))
        self.assertEqual(count, 1)
        np.testing.assert_allclose(geometry.bounds, [-.5, -3.5, 1.5, 4.5])
        # A ring of seeded bin centers makes a hole, retained in the overlay source.
        ring = [f'b{x}_{y}' for y in range(3) for x in range(3) if x in (0, 2) or y in (0, 2)]
        ring_group = create_seed_group(layer, ring, 'Ring')
        shape, n = manager.compute_group_radius(ring_group, .6, (1., 1.))
        self.assertEqual(n, 8)
        self.assertEqual(shape.geom_type, 'Polygon')
        self.assertGreater(len(shape.interiors), 0)
        manager.set_draft(shape, 'radius')
        self.assertGreater(len(manager.source.data['xs'][0][0]), 1)
        document = Document()
        document.add_root(manager.component.figure)
        document.to_json()

    def test_latest_preview_wins_and_cancel_discards_work(self):
        layer = make_layer()
        manager = layer.roi_manager
        group = create_seed_group(layer, ['b0_0', 'b1_0'])
        manager.method.value = 'radius'
        manager.seed_source.value = 'group'
        manager.seed_group.value = group
        manager.radius.value = .2
        real = manager.compute_group_radius
        running, release = threading.Event(), threading.Event()
        radii = []
        def slow(group, radius, *args):
            radii.append(radius)
            if len(radii) == 1:
                running.set()
                release.wait(10)
            return real(group, radius, *args)
        manager.compute_group_radius = slow
        manager.preview_group_radius(immediate=True)
        self.assertTrue(running.wait(5))
        manager.radius.value = .8
        self.assertTrue(manager.save_button.disabled)
        release.set()
        wait_for(lambda: not manager.logger.running and manager._pending_preview is None)
        self.assertEqual(radii, [.2, .8])
        expected, _ = real(group, .8, (1., 1.))
        self.assertLess(manager.draft.symmetric_difference(expected).area, 1e-10)
        manager.radius.value = .4
        manager.cancel()
        self.assertIsNone(manager._preview_timeout)
        self.assertIsNone(manager._pending_preview)
        self.assertIsNone(manager.draft)
        self.assertTrue(manager.save_button.disabled)

    def test_group_removal_invalid_inputs_and_manual_seed_mode(self):
        layer = make_layer()
        manager = layer.roi_manager
        group = create_seed_group(layer, ['b0_0'])
        manager.method.value = 'radius'
        self.assertTrue(manager.seed_source.visible)
        manager.seed_source.value = 'group'
        manager.seed_group.value = group
        manager.radius.value = 0
        self.assertIn('positive', manager.logger.div.text)
        self.assertIsNone(manager.draft)
        manager.radius.value = .2
        fire_preview(manager)
        layer.groups_manager.delete_group(group, manager.doc)
        self.assertEqual(manager.seed_group.value, '')
        self.assertIsNone(manager.draft)
        self.assertTrue(manager.save_button.disabled)
        manager.seed_source.value = 'map'
        self.assertEqual(manager.draw_button.label, 'Draw')
        self.assertFalse(manager.seed_group.visible)
        manager.begin()
        self.assertIs(manager.drawing.owner, manager)
        manager.cancel()
        for kind in ('spots', 'other_cells'):
            other = make_layer(kind).roi_manager
            other.method.value = 'radius'
            self.assertFalse(other.seed_source.visible)
            self.assertFalse(other._group_mode())


if __name__ == '__main__':
    unittest.main()
