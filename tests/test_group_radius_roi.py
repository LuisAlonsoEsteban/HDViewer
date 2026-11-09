"""ROIs seeded by the bins of a group: a radius around each bin, optionally minus the group."""
import threading
import unittest
from unittest.mock import patch

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import union_all
from shapely.geometry import MultiPolygon, Polygon, box

from hd_viewer import settings
from hd_viewer.core.roi import ROIGeometry
from hd_viewer.startup.sample import Bins
from tests.test_roi_manager import Api, make_layer

B1_1 = 5            # native index of bin b1_1


def create_seed_group(view, names, name='Seeds'):
    return view.groups.add_fixed_group(name, view.sample.adata.obs_names.isin(names))


def group_shape(group, radius, exclude_group=False, **extra):
    return dict(kind='group', group=group, radius=radius, exclude_group=exclude_group, **extra)


def outline(xs, ys):
    """Geometry of a multi_polygons outline (as sent to the browser)."""
    polygons = [Polygon(list(zip(rings_x[0], rings_y[0])),
                        [list(zip(hx, hy)) for hx, hy in zip(rings_x[1:], rings_y[1:])])
                for rings_x, rings_y in zip(xs, ys)]
    return polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)


class GroupRadiusGeometryTests(unittest.TestCase):
    def test_complex_overlay_is_bounded_and_keeps_exact_geometry_and_membership(self):
        api = Api(self)
        rois = api.view.rois
        group = create_seed_group(api.view, ['b1_1'])
        with patch.object(settings, 'ROI_OVERLAY_MAX_VERTICES', 32), \
             patch.object(settings, 'ROI_OVERLAY_MAX_SIZE', 96):
            geometry, count = rois.compute_group_radius(group, 1.2, (1., 1.), exclude_group=True)
            self.assertEqual(count, 1)
            exact = rois.radius_geometry(1.5, 1.5, 1.2).difference(box(1, 1, 2, 2))
            self.assertLess(geometry.symmetric_difference(exact).area, 1e-10)
            overlay = rois._prepare_overlay(geometry)
            self.assertIsNone(overlay.xs)
            image = overlay.image('#ff0000')
            self.assertEqual(image.dtype, np.uint32)
            self.assertLessEqual(image.nbytes, 96 * 96 * 4)
            rgba = image.view(np.uint8).reshape((*image.shape, 4))
            x0, y0, x1, y1 = overlay.bounds

            def alpha_at(x, y):
                ix = int((x - x0) / (x1 - x0) * image.shape[1])
                iy = int((y - y0) / (y1 - y0) * image.shape[0])
                return rgba[iy, ix, 3]
            self.assertEqual(alpha_at(1.5, 1.5), 0)  # Source-group hole stays transparent.
            self.assertGreater(alpha_at(1.5, 2.2), 0)

            # Through the API: the outline sent is bounded, the saved shape and membership exact
            preview = api.ok(api.post('rois/preview', json=group_shape(group, 1.2, True)))
            sent = sum(len(ring) for rings in preview['xs'] for ring in rings)
            self.assertLessEqual(sent, 32)
            saved = api.roi(**group_shape(group, 1.2, True))
        roi = rois.rois[0]
        self.assertLess(roi.geometry.symmetric_difference(exact).area, 1e-10)
        expected_membership = rois.compute_membership(exact)
        expected_membership[B1_1] = False
        np.testing.assert_array_equal(roi.mask, expected_membership)
        self.assertEqual(saved['count'], int(expected_membership.sum()))

    def test_raster_dimensions_cannot_exceed_budget_from_float_roundoff(self):
        with patch.object(settings, 'ROI_OVERLAY_MAX_VERTICES', 1):
            overlay = ROIGeometry._prepare_overlay(box(-.9, -.9, 189.9, 189.9))
        self.assertLessEqual(max(overlay.raster.shape), settings.ROI_OVERLAY_MAX_SIZE)

    def test_display_preparation_and_group_radius_can_be_cancelled(self):
        view = make_layer()
        group = create_seed_group(view, ['b1_1'])
        geometry, _ = view.rois.compute_group_radius(group, 1.2, (1., 1.))
        with self.assertRaises(InterruptedError):
            ROIGeometry._prepare_overlay(geometry, lambda: True)
        for exclude in (False, True):
            with self.subTest(exclude=exclude), self.assertRaises(InterruptedError):
                view.rois.compute_group_radius(group, 1.2, (1., 1.), cancelled=lambda: True, exclude_group=exclude)

    def test_source_group_union_is_reused_for_radius_changes_and_invalidated_for_membership(self):
        view = make_layer()
        rois = view.rois
        group = create_seed_group(view, ['b1_1'])
        rois.compute_group_radius(group, 1.2, (1., 1.), exclude_group=True)
        cached = rois._group_union_cache[3]
        rois.compute_group_radius(group, 1.5, (1., 1.), exclude_group=True)
        self.assertIs(rois._group_union_cache[3], cached)
        view.sample.adata.obs[f'grp_{group}'] = view.sample.adata.obs_names == 'b3_3'
        rois.compute_group_radius(group, 1.5, (1., 1.), exclude_group=True)
        self.assertIsNot(rois._group_union_cache[3], cached)
        self.assertTrue(rois._group_union_cache[3].equals(box(3, 3, 4, 4)))

    def test_large_group_uses_every_seed_across_union_batches(self):
        n = 23
        names = [f'b{x}_{y}' for y in range(n) for x in range(n)]
        gdf = gpd.GeoDataFrame({'index': names}, index=names,
                               geometry=[box(x, y, x + 1, y + 1) for y in range(n) for x in range(n)])
        data = ad.AnnData(np.ones((len(names), 1)), obs=pd.DataFrame(index=names),
                          var=pd.DataFrame(index=['g0']))
        data.obs['grp_All'] = True
        rois = ROIGeometry(Bins(gdf, data, pixels_per_bin=1, bin_size=2))
        shape, count = rois.compute_group_radius('All', .8, (1., 1.))
        self.assertEqual(count, 529)
        self.assertTrue(shape.is_valid)
        np.testing.assert_allclose(shape.bounds, [-.3, -.3, 23.3, 23.3])
        self.assertTrue(rois.compute_membership(shape).all())

    def test_minus_group_with_small_radius_explains_empty_result(self):
        api = Api(self)
        group = create_seed_group(api.view, ['b0_0'])
        with self.assertRaisesRegex(ValueError, 'Increase the radius'):
            api.view.rois.compute_group_radius(group, .2, (1., 1.), exclude_group=True)
        response = api.post('rois/preview', json=group_shape(group, .2, True))
        self.assertEqual(response.status_code, 422)
        self.assertIn('leaves no ROI', response.json()['detail'])
        self.assertEqual(api.post('rois', json=group_shape(group, .2, True)).status_code, 422)
        self.assertEqual(api.ok(api.get('rois'))['rois'], [])
        self.assertEqual(api.ok(api.post('rois/preview', json=group_shape(group, .2)))['seeds'], 1)

    def test_subtraction_uses_native_shapes_after_radius_calibration(self):
        view = make_layer()
        group = create_seed_group(view, ['b0_0', 'b3_3'])
        full, _ = view.rois.compute_group_radius(group, 2., (2., .5))
        result, count = view.rois.compute_group_radius(group, 2., (2., .5), exclude_group=True)
        expected = full.difference(union_all([box(0, 0, 1, 1), box(3, 3, 4, 4)]))
        self.assertEqual(count, 2)
        self.assertTrue(result.is_valid)
        self.assertEqual(result.geom_type, 'MultiPolygon')
        self.assertLess(result.symmetric_difference(expected).area, 1e-10)

    def test_physical_scale_and_holes_are_preserved(self):
        view = make_layer()
        rois = view.rois
        group = create_seed_group(view, ['b0_0'])
        geometry, count = rois.compute_group_radius(group, 2., (2., .5))
        self.assertEqual(count, 1)
        np.testing.assert_allclose(geometry.bounds, [-.5, -3.5, 1.5, 4.5])
        # A ring of seeded bin centers makes a hole, retained in the outline.
        ring = [f'b{x}_{y}' for y in range(3) for x in range(3) if x in (0, 2) or y in (0, 2)]
        ring_group = create_seed_group(view, ring, 'Ring')
        shape, n = rois.compute_group_radius(ring_group, .6, (1., 1.))
        self.assertEqual(n, 8)
        self.assertEqual(shape.geom_type, 'Polygon')
        self.assertGreater(len(shape.interiors), 0)
        overlay = rois._prepare_overlay(shape)
        self.assertGreater(len(overlay.xs[0]), 1)

    def test_categorical_group_conditions_work(self):
        view = make_layer()
        view.sample.adata.obs['cell_type'] = ['Tumor', 'Tumor', 'Normal', 'Normal'] * 4
        mask = view.groups._build_mask([('cell_type', '==', 'Tumor', 'AND')], view.sample.adata)
        np.testing.assert_array_equal(mask.to_numpy(), np.tile([True, True, False, False], 4))


class GroupRadiusApiTests(unittest.TestCase):
    def test_all_group_bins_create_one_roi_even_outside_viewport_and_at_aggregate_level(self):
        api = Api(self)
        view = api.view
        group = create_seed_group(view, ['b0_0', 'b3_3'])
        # Hide the group and load an aggregated viewport away from one seed.
        self.assertEqual(api.patch(f'groups/{group}', json=dict(visible=False)).status_code, 200)
        self.assertEqual(view.groups._active_groups, set())
        view.rows.set_fetch_context(bbox=(0, 1, 0, 1), resolution_level='agg:2')
        view.rows.fetch_data()
        self.assertIn(group, api.ok(api.get('rois'))['seed_groups'])
        preview = api.ok(api.post('rois/preview', json=group_shape(group, .2)))
        self.assertEqual(preview['seeds'], 2)
        self.assertEqual(len(preview['xs']), 2)                 # One ROI, with two pieces.
        drawn = outline(preview['xs'], preview['ys'])
        expected = union_all([ROIGeometry.radius_geometry(.5, .5, .2), ROIGeometry.radius_geometry(3.5, 3.5, .2)])
        self.assertLess(drawn.symmetric_difference(expected).area, 1e-10)
        self.assertEqual(api.ok(api.get('rois'))['rois'], [])   # a preview saves nothing
        self.assertEqual(view.selected, set())

        roi = api.roi(**group_shape(group, .2, name='Group neighbourhood'))
        self.assertEqual(len(api.ok(api.get('rois'))['rois']), 1)
        self.assertEqual((roi['name'], roi['seed_count'], roi['seed_group'], roi['count']),
                         ('Group neighbourhood', 2, group, 2))
        self.assertFalse(roi['group_excluded'])
        saved = view.rois.rois[0]
        self.assertEqual(saved.geometry.geom_type, 'MultiPolygon')
        self.assertLess(saved.geometry.symmetric_difference(expected).area, 1e-10)
        created = api.roi_group(roi)
        np.testing.assert_array_equal(view.sample.adata.obs[f'grp_{created}'], saved.mask)
        api.ok(api.delete(f"rois/{roi['id']}"))
        self.assertIn(group, view.groups._groups)
        self.assertNotIn(created, view.groups._groups)

    def test_minus_group_preview_geometry_membership_and_saved_snapshot(self):
        api = Api(self)
        view, rois = api.view, api.view.rois
        group = create_seed_group(view, ['b1_1'])
        full_preview = api.ok(api.post('rois/preview', json=group_shape(group, 1.2)))
        full = outline(full_preview['xs'], full_preview['ys'])
        full_members = rois.compute_membership(full)
        self.assertTrue(full_members[B1_1])
        minus_preview = api.ok(api.post('rois/preview', json=group_shape(group, 1.2, True)))
        minus = outline(minus_preview['xs'], minus_preview['ys'])
        expected = full.difference(box(1, 1, 2, 2))
        self.assertLess(minus.symmetric_difference(expected).area, 1e-10)
        self.assertEqual(len(minus.interiors), 1)
        self.assertGreater(len(minus_preview['xs'][0]), 1)                 # the hole is drawn
        # Normal intersection still includes the source bin along the hole edge.
        self.assertTrue(rois.compute_membership(minus)[B1_1])

        roi = api.roi(**group_shape(group, 1.2, True))
        self.assertTrue(roi['group_excluded'])
        saved = rois.rois[0]
        expected_members = full_members.copy()
        expected_members[B1_1] = False
        np.testing.assert_array_equal(saved.mask, expected_members)
        self.assertEqual(saved.seed_group, group)
        api.ok(api.post(f"rois/{roi['id']}/select"))
        self.assertNotIn('b1_1', view.selected)
        self.assertEqual(len(view.selected), int(expected_members.sum()))
        created = api.roi_group(roi)
        np.testing.assert_array_equal(view.sample.adata.obs[f'grp_{created}'], expected_members)
        # Later source edits and previews do not modify saved ROIs.
        view.sample.adata.obs[f'grp_{group}'] = view.sample.adata.obs_names == 'b3_3'
        moved = api.ok(api.post('rois/preview', json=group_shape(group, 1.2)))
        np.testing.assert_allclose(outline(moved['xs'], moved['ys']).centroid.coords[0], [3.5, 3.5], atol=1e-6)
        np.testing.assert_array_equal(saved.mask, expected_members)
        self.assertLess(saved.geometry.symmetric_difference(expected).area, 1e-10)

    def test_radius_changes_rebuild_preview_without_modifying_saved_rois(self):
        api = Api(self)
        view = api.view
        group = create_seed_group(view, ['b0_0', 'b1_0'])
        self.assertEqual(len(api.ok(api.post('rois/preview', json=group_shape(group, .2)))['xs']), 2)
        api.roi(**group_shape(group, .2))
        saved = view.rois.rois[0]
        self.assertEqual(saved.geometry.geom_type, 'MultiPolygon')
        saved_area = saved.geometry.area
        bigger = api.ok(api.post('rois/preview', json=group_shape(group, 1.)))
        self.assertEqual(len(bigger['xs']), 1)                     # one polygon
        self.assertGreater(outline(bigger['xs'], bigger['ys']).area, saved_area)
        self.assertEqual(saved.geometry.area, saved_area)
        self.assertEqual(len(api.ok(api.get('rois'))['rois']), 1)
        # Membership changes of the seed group update the preview too.
        view.sample.adata.obs[f'grp_{group}'] = view.sample.adata.obs_names == 'b3_3'
        moved = api.ok(api.post('rois/preview', json=group_shape(group, 1.)))
        self.assertEqual(moved['seeds'], 1)
        np.testing.assert_allclose(outline(moved['xs'], moved['ys']).centroid.coords[0], [3.5, 3.5], atol=1e-6)

    def test_latest_preview_wins(self):
        for changed in ('exclude_group', 'radius'):
            with self.subTest(changed=changed):
                api = Api(self)
                rois = api.view.rois
                group = create_seed_group(api.view, ['b0_0', 'b1_0'])
                first = group_shape(group, 1.2)
                second = dict(first, exclude_group=True) if changed == 'exclude_group' else dict(first, radius=.8)
                real = rois.compute_group_radius
                entered, release = threading.Event(), threading.Event()
                calls = []

                def slow(group, radius, *args, **kwargs):
                    calls.append((radius, kwargs['exclude_group']))
                    if len(calls) == 1:
                        entered.set()
                        release.wait(10)
                    return real(group, radius, *args, **kwargs)
                rois.compute_group_radius = slow
                responses = {}
                worker = threading.Thread(target=lambda: responses.setdefault(
                    'first', api.post('rois/preview', json=first)))
                worker.start()
                self.assertTrue(entered.wait(10))
                latest = api.ok(api.post('rois/preview', json=second))
                release.set()
                worker.join(30)
                self.assertEqual(responses['first'].status_code, 409)   # obsolete: discarded
                self.assertEqual(calls, [(1.2, False), (second['radius'], second['exclude_group'])])
                expected, _ = real(group, second['radius'], (1., 1.), exclude_group=second['exclude_group'])
                self.assertLess(outline(latest['xs'], latest['ys']).symmetric_difference(expected).area, 1e-10)
                # Saving the latest preview reuses its geometry
                api.roi(**second)
                self.assertEqual(len(calls), 2)
                self.assertLess(rois.rois[0].geometry.symmetric_difference(expected).area, 1e-10)
                self.assertEqual(rois.rois[0].excluded_seed_mask is not None, second['exclude_group'])

    def test_group_removal_invalid_inputs_and_other_layers(self):
        api = Api(self)
        group = create_seed_group(api.view, ['b0_0'])
        response = api.post('rois/preview', json=group_shape(group, 0))
        self.assertEqual(response.status_code, 422)
        self.assertIn('positive', response.json()['detail'])
        self.assertEqual(api.post('rois/preview', json=dict(kind='group', group=group)).status_code, 422)
        api.ok(api.post('rois/preview', json=group_shape(group, .2)))
        api.ok(api.delete(f'groups/{group}'))
        self.assertEqual(api.ok(api.get('rois'))['seed_groups'], [])
        response = api.post('rois/preview', json=group_shape(group, .2))
        self.assertEqual(response.status_code, 422)
        self.assertIn('no longer exists', response.json()['detail'])
        # Seeds from a map point still work
        self.assertEqual(api.roi(kind='point', centre=[.5, .5], radius=.2)['count'], 1)
        for kind in ('spots', 'cells'):
            with self.subTest(kind=kind):
                other = Api(self, kind)
                name = other.view.groups.add_fixed_group('Seeds', np.arange(16) == 0)
                self.assertEqual(other.ok(other.get('rois'))['seed_groups'], [])
                response = other.post('rois/preview', json=group_shape(name, .2))
                self.assertEqual(response.status_code, 422)
                self.assertIn('bins layers', response.json()['detail'])


if __name__ == '__main__':
    unittest.main()
