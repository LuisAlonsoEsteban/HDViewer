"""Regions of interest: geometry, exact membership, groups, DE and selection.

A 4 x 4 grid of unit squares (bin b{x}_{y} at [x, x+1] x [y, y+1], native index
4y + x) plus one small not-sequenced shape inside b1_1.
"""
import threading
import time
import unittest
from types import SimpleNamespace

import numpy as np
from shapely.geometry import box

from hd_viewer.analysis.de import DEAnalysis
from hd_viewer.core.encoding import StaticParameter
from hd_viewer.core.groups import GROUPS_COL
from hd_viewer.core.roi import ROI, ROIGeometry
from hd_viewer.core.selection import aggregation_names
from tests.fixtures import make_sample, make_view, show
from tests.test_clustering_markers import api_client

KINDS = ('bins', 'spots', 'cells')


def make_layer(kind='bins'):
    """The API view of the 4 x 4 grid, with every row loaded."""
    view = make_view(kind=kind, n=4, not_sequenced=1)
    show(view, -1, 5, -1, 5)
    return view


def saved(view, name='Region', geometry=None):
    geometry = box(.2, .2, 1.8, 1.8) if geometry is None else geometry
    roi = ROI(name, geometry, 'box', '#ff0000')
    roi.mask = view.rois.membership(roi)
    view.rois.rois.append(roi)
    return roi


def wait_for(predicate, message='Timed out'):
    deadline = time.monotonic() + 60
    while not predicate() and time.monotonic() < deadline:
        time.sleep(.01)
    if not predicate():
        raise AssertionError(message)


class Api:
    """The grid behind the HTTP API (layer id = kind), every row loaded."""
    def __init__(self, test, kind='bins'):
        self.test, self.kind = test, kind
        self.client, self.session = api_client(test, make_sample(kind, n=4, not_sequenced=1),
                                               kind=kind, layer_id=kind)
        self.view = self.session.views[kind]
        self.post('viewport', json=dict(seq=1, x0=-1, x1=5, y0=-1, y1=5, width=900, height=900))

    def url(self, path):
        return f'/api/layers/{self.kind}/{path}'

    def get(self, path, **kwargs):
        return self.client.get(self.url(path), **kwargs)

    def post(self, path, **kwargs):
        return self.client.post(self.url(path), **kwargs)

    def patch(self, path, **kwargs):
        return self.client.patch(self.url(path), **kwargs)

    def delete(self, path, **kwargs):
        return self.client.delete(self.url(path), **kwargs)

    def ok(self, response):
        self.test.assertEqual(response.status_code, 200, response.text)
        return response.json()

    def roi(self, **shape):
        return self.ok(self.post('rois', json=shape))['rois'][-1]

    def roi_group(self, roi):
        return self.ok(self.post(f"rois/{roi['id']}/group"))['group']


class ROIGeometryTests(unittest.TestCase):
    def test_membership_uses_intersection_and_full_native_data(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                view = make_layer(kind)
                rois = view.rois
                # Nine native shapes touch or intersect the box; an unsequenced shape is excluded.
                mask = rois.compute_membership(box(1, 1, 2, 2))
                self.assertEqual(mask.shape, (16,))
                self.assertEqual(int(mask.sum()), 9)
                show(view, 0, 1, 0, 1)                          # rows far from the ROI not needed
                far = saved(view, 'Far', box(2.1, 2.1, 3.9, 3.9))
                self.assertEqual(far.count, 4)
                self.assertEqual(rois.compute_membership(box(10, 10, 11, 11)).sum(), 0)
                lasso = rois.lasso_geometry([.2, 1.8, 1.8, .2], [.2, .2, 1.8, 1.8])
                np.testing.assert_array_equal(rois.compute_membership(lasso), saved(view).mask)

    def test_geometry_validation_and_unequal_radius_scales(self):
        invalid = [lambda: ROIGeometry.box_geometry(0, 0, 0, 1),
                   lambda: ROIGeometry.box_geometry(0, 0, np.inf, 1),
                   lambda: ROIGeometry.lasso_geometry([0, 1], [0, 1]),
                   lambda: ROIGeometry.lasso_geometry([0, 1, 0, 1], [0, 1, 1, 0]),
                   lambda: ROIGeometry.lasso_geometry([0, 1, np.nan], [0, 1, 2]),
                   lambda: ROIGeometry.radius_geometry(0, 0, 0),
                   lambda: ROIGeometry.radius_geometry(0, 0, -1)]
        for draw in invalid:
            with self.assertRaises(ValueError):
                draw()
        shape = ROIGeometry.radius_geometry(5, 6, 2, 2, .5)
        np.testing.assert_allclose(shape.bounds, (4, 2, 6, 10))
        np.testing.assert_allclose(ROIGeometry.box_geometry(2, 3, 0, 1).bounds, (0, 1, 2, 3))

    def test_focus_visibility_names_and_empty_region(self):
        api = Api(self)
        roi = api.roi(kind='box', name='Region', box=[.2, .2, 1.8, 1.8])
        self.assertEqual(roi['count'], 4)
        focus = api.ok(api.post(f"rois/{roi['id']}/focus"))
        np.testing.assert_allclose([focus['x0'], focus['y0'], focus['x1'], focus['y1']], [.12, .12, 1.88, 1.88])
        listed = api.ok(api.patch(f"rois/{roi['id']}", json=dict(visible=False)))['rois']
        self.assertFalse(listed[0]['visible'])
        listed = api.ok(api.patch(f"rois/{roi['id']}", json=dict(visible=True)))['rois']
        self.assertTrue(listed[0]['visible'])
        duplicate = api.post('rois', json=dict(kind='box', name='Region', box=[.1, .1, .9, .9]))
        self.assertEqual(duplicate.status_code, 422)
        self.assertIn('already exists', duplicate.json()['detail'])

        empty = api.roi(kind='box', name='Empty', box=[10, 10, 11, 11])
        self.assertEqual(empty['count'], 0)
        self.assertEqual(api.post(f"rois/{empty['id']}/group").status_code, 422)
        self.assertEqual(api.ok(api.post(f"rois/{empty['id']}/select"))['selected'], 0)
        self.assertEqual(api.post(f"rois/{empty['id']}/focus").status_code, 200)
        remaining = api.ok(api.delete(f"rois/{empty['id']}"))['rois']
        self.assertEqual([r['name'] for r in remaining], ['Region'])
        self.assertEqual(api.delete(f"rois/{empty['id']}").status_code, 404)

    def test_seed_radius_uses_calibration_or_full_resolution_pixels(self):
        view = make_layer()
        rois = view.rois
        rois.coordinate_scale = SimpleNamespace(data=dict(px_x=[2.], px_y=[1.], um_per_px=[.5]))
        factors, unit = rois.radius_scale()
        self.assertEqual(unit, 'µm')
        np.testing.assert_allclose(rois.seed_radius_geometry(2, 2, 2., *factors).bounds, (0, -2, 4, 6))
        rois.coordinate_scale = SimpleNamespace(data=dict(px_x=[2.], px_y=[4.], um_per_px=[None]))
        factors, unit = rois.radius_scale()
        self.assertEqual(unit, 'full-res px')
        np.testing.assert_allclose(rois.seed_radius_geometry(2, 2, 2., *factors).bounds, (1, 1.5, 3, 2.5))
        with self.assertRaisesRegex(ValueError, 'radius'):
            rois.seed_radius_geometry(2, 2, -1., *factors)

        api = Api(self)
        self.assertEqual(api.ok(api.get('rois'))['unit'], 'full-res px')
        response = api.post('rois', json=dict(kind='point', centre=[2, 2], radius=-1))
        self.assertEqual(response.status_code, 422)
        self.assertIn('radius', response.json()['detail'])
        self.assertEqual(api.ok(api.get('rois'))['rois'], [])


class ROIGroupTests(unittest.TestCase):
    def test_de_in_flight_is_discarded_when_linked_roi_is_deleted(self):
        api = Api(self)
        roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
        name = api.roi_group(roi)
        analysis = DEAnalysis(api.view.sample)
        real = analysis.run
        running, release = threading.Event(), threading.Event()

        def blocked(*args, **kwargs):
            running.set()
            release.wait(10)
            return real(*args, **kwargs)
        analysis.run = blocked
        api.view.de['analysis'] = analysis
        self.assertEqual(api.post('de', json=dict(groups=[name])).status_code, 202)
        self.assertTrue(running.wait(5))
        api.ok(api.delete(f"rois/{roi['id']}"))
        release.set()
        wait_for(lambda: api.view.job['state'] != 'running')
        self.assertIn('discarded', api.view.job['message'])
        self.assertIsNone(api.view.de['header'])
        self.assertEqual(api.get('de/table').status_code, 404)

    def test_fixed_groups_survive_conditions_and_delete_with_roi(self):
        api = Api(self)
        obs = api.view.sample.adata.obs
        roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
        group_name = api.roi_group(roi)
        mask = api.view.rois.rois[0].mask
        np.testing.assert_array_equal(obs[f'grp_{group_name}'], mask)
        api.ok(api.post('groups', json=dict(name='Condition', conditions=[dict(column='g0', op='>', value='15')])))
        np.testing.assert_array_equal(obs[f'grp_{group_name}'], mask)
        self.assertIn(group_name, [g['name'] for g in api.ok(api.get('de/options'))['groups']])
        api.view.de.update(header='DE vs Rest', merged=[], compared=('groups', [group_name]))
        removed = api.ok(api.delete(f"rois/{roi['id']}"))
        self.assertTrue(removed['recoloured'])
        self.assertEqual([g['name'] for g in removed['groups']], ['Condition'])
        self.assertNotIn(group_name, api.view.groups._groups)
        self.assertNotIn(f'grp_{group_name}', obs)
        self.assertIsNone(api.view.de['header'])
        self.assertIn('Condition', api.view.groups._groups)

    def test_manual_group_delete_recreation_and_unique_names(self):
        api = Api(self)
        groups = api.view.groups
        roi = api.roi(kind='box', name='Region', box=[.2, .2, 1.8, 1.8])
        groups.add_fixed_group('ROI: Region', np.arange(16) < 4)
        linked = api.roi_group(roi)
        self.assertEqual(linked, 'ROI: Region (2)')
        self.assertEqual(api.roi_group(roi), linked)
        api.ok(api.delete(f'groups/{linked}'))
        self.assertIsNone(groups.group_for_roi(roi['id']))
        self.assertIsNone(api.ok(api.get('rois'))['rois'][0]['group'])
        self.assertEqual(api.roi_group(roi), linked)
        api.ok(api.delete(f"rois/{roi['id']}"))
        self.assertIn('ROI: Region', groups._groups)
        api.ok(api.delete('groups/ROI: Region'))
        self.assertEqual(groups._groups, {})
        self.assertNotIn(GROUPS_COL, api.view.sample.adata.obs)
        self.assertIsInstance(api.view.rows.get_param('color'), StaticParameter)
        self.assertNotEqual(api.view.colour_by['kind'], 'groups')

    def test_roi_group_runs_de(self):
        api = Api(self)
        roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
        name = api.roi_group(roi)
        self.assertEqual(api.post('de', json=dict(groups=[name])).status_code, 202)
        wait_for(lambda: api.view.job['state'] != 'running')
        self.assertEqual(api.view.job['state'], 'done', api.view.job['message'])
        self.assertEqual(api.view.de['names'], [name])
        api.ok(api.delete(f"rois/{roi['id']}"))
        self.assertIsNone(api.view.de['header'])


class ROISelectionTests(unittest.TestCase):
    def test_native_selection_projects_to_aggregated_map_and_embedding(self):
        api = Api(self)
        view, sample = api.view, api.view.sample
        roi = api.roi(kind='box', box=[.1, .1, .9, .9])
        rows = view.rows
        rows.set_fetch_context(bbox=(-1, 5, -1, 5), resolution_level='agg:2')
        rows.fetch_data()
        names = list(rows.loaded_obs_names)
        view.embedding.update(method='pca', names=np.asarray(sample.adata.obs_names.astype(str)))

        selected = api.ok(api.post(f"rois/{roi['id']}/select"))
        self.assertEqual(view.selected, {'b0_0'})
        expected = names.index('summary:0:0:2')
        self.assertEqual(selected['rows'], [expected])
        self.assertEqual(selected['row_version'], rows.row_version)
        self.assertEqual(api.ok(api.get('embedding/selection'))['indices'], [0])
        # Selecting an aggregate expands to its native members; native embeddings see all of them.
        self.assertEqual(api.ok(api.post('selection', json=dict(row_version=rows.row_version,
                                                                  rows=[expected])))['selected'], 4)
        self.assertEqual(view.selected, {'b0_0', 'b1_0', 'b0_1', 'b1_1'})
        self.assertEqual(api.ok(api.get('selection'))['rows'], [expected])
        self.assertEqual(api.ok(api.get('embedding/selection'))['indices'], [0, 1, 4, 5])
        # Stale row positions are refused
        self.assertEqual(api.post('selection', json=dict(row_version=rows.row_version - 1, rows=[0])).status_code, 409)
        # Projection did not need expression aggregation.
        self.assertEqual(sample._aggregated_samples, {})
        self.assertEqual(aggregation_names(sample, 2).loc['b0_0'], 'summary:0:0:2')
        rows.set_fetch_context(bbox=(-1, 5, -1, 5), resolution_level=None)
        rows.fetch_data()
        positions = api.ok(api.get('selection'))['rows']
        self.assertEqual(sorted(rows.loaded_obs_names[i] for i in positions), ['b0_0', 'b0_1', 'b1_0', 'b1_1'])
        # Points chosen on the embedding become the map's selection
        api.ok(api.post('embedding/selection', json=dict(indices=[15])))
        self.assertEqual(view.selected, {'b3_3'})


if __name__ == '__main__':
    unittest.main()
