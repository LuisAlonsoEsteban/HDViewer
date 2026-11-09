"""Groups on the map (visibility) and the layer's selection when ROIs and groups change."""
import unittest

import numpy as np

from hd_viewer.core.groups import GROUPS_COL
from tests.fixtures import loaded
from tests.test_roi_manager import KINDS, Api, make_layer


class GroupSelectionTests(unittest.TestCase):
    def test_group_visibility_shows_only_the_shown_groups(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                view = make_layer(kind)
                groups = view.groups
                obs_names = view.sample.adata.obs_names
                first, last = obs_names[0], obs_names[-1]
                name = groups.add_fixed_group('A', obs_names == first)
                groups.refresh_active_labels()
                view.colour_by_categories(GROUPS_COL, groups._color_map, 'groups')

                def visible():
                    mask = view.visibility()
                    self.assertIsNotNone(mask)
                    return dict(zip(loaded(view), mask.astype(bool)))

                shown = visible()
                self.assertTrue(shown[first])
                self.assertFalse(shown[last])
                self.assertEqual(sum(shown.values()), 1)            # only the shown group's member
                view.rows.clear_sources()                           # reloaded rows keep the rule
                view.rows.fetch_data()
                self.assertEqual(sum(visible().values()), 1)
                groups._active_groups.discard(name)
                groups.refresh_active_labels()
                self.assertEqual(sum(visible().values()), 0)

    def test_clear_selection(self):
        api = Api(self)
        roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
        self.assertEqual(api.ok(api.post(f"rois/{roi['id']}/select"))['selected'], 4)
        cleared = api.ok(api.post('selection', json=dict(row_version=api.view.rows.row_version, rows=[])))
        self.assertEqual(cleared['selected'], 0)
        self.assertEqual(api.view.selected, set())
        self.assertEqual(api.ok(api.get('selection'))['rows'], [])

    def test_unselect_and_remove_preserve_unrelated_selection(self):
        for action in ('unselect', 'remove_roi', 'delete_group'):
            with self.subTest(action=action):
                api = Api(self)
                view = api.view
                roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
                name = api.roi_group(roi)
                members = set(view.sample.adata.obs_names[view.rois.rois[0].mask])
                view.selected = members | {'b3_3'}
                if action == 'delete_group':
                    api.ok(api.delete(f'groups/{name}'))
                elif action == 'remove_roi':
                    api.ok(api.delete(f"rois/{roi['id']}"))
                else:
                    api.ok(api.post(f"rois/{roi['id']}/unselect"))
                self.assertEqual(view.selected, {'b3_3'})
                displayed = loaded(view)
                self.assertEqual({displayed[i] for i in api.ok(api.get('selection'))['rows']}, {'b3_3'})

    def test_deletion_clears_embedding_selection_too(self):
        api = Api(self)
        view = api.view
        roi = api.roi(kind='box', box=[.2, .2, 1.8, 1.8])
        group = api.roi_group(roi)
        view.embedding.update(method='pca', names=np.asarray(view.sample.adata.obs_names.astype(str)))
        api.ok(api.post(f"rois/{roi['id']}/select"))
        self.assertEqual(api.ok(api.get('embedding/selection'))['indices'], [0, 1, 4, 5])
        api.ok(api.delete(f'groups/{group}'))
        self.assertEqual(api.ok(api.get('embedding/selection'))['indices'], [])
        self.assertEqual(api.ok(api.get('selection'))['rows'], [])

    def test_condition_creation_failure_creates_no_group(self):
        api = Api(self)
        response = api.post('groups', json=dict(conditions=[dict(column='g0', op='>', value='invalid')]))
        self.assertEqual(response.status_code, 422)
        self.assertIn('numeric', response.json()['detail'])
        self.assertEqual(api.ok(api.get('groups'))['groups'], [])
        self.assertNotIn('grp_Group A', api.view.sample.adata.obs)
        for bad in (dict(column='nope', op='==', value='1'), dict(column='g0', op='!=', value='1')):
            with self.subTest(condition=bad):
                self.assertEqual(api.post('groups/check', json=bad).status_code, 422)
        self.assertEqual(api.post('groups', json=dict(conditions=[])).status_code, 422)

    def test_condition_creation_and_deletion(self):
        api = Api(self)
        hint = api.ok(api.post('groups/check', json=dict(column='g0', op='>', value='15')))['hint']
        self.assertIn('max=31.00', hint)
        created = api.ok(api.post('groups', json=dict(conditions=[dict(column='g0', op='>', value='15')])))
        self.assertTrue(created['shown'])
        self.assertEqual([(g['name'], g['size']) for g in created['groups']], [('Group A', 8)])
        self.assertEqual(int(api.view.sample.adata.obs['grp_Group A'].sum()), 8)
        self.assertEqual(api.post('groups', json=dict(name='Group A', conditions=[
            dict(column='g0', op='>', value='1')])).status_code, 422)          # names are unique
        api.view.selected = {'b3_3'}
        api.ok(api.delete('groups/Group A'))
        self.assertEqual(api.view.selected, set())
        self.assertEqual(api.ok(api.get('groups'))['groups'], [])


if __name__ == '__main__':
    unittest.main()
