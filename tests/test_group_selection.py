import threading
import unittest
from unittest.mock import patch

from bokeh.document import Document
from bokeh.server.callbacks import NextTickCallback
from bokeh.events import ButtonClick
from bokeh.models import Button

from hd_viewer.startup.sample import Bins  # Initialize startup before importing board components.
from hd_viewer.board.layout.embedding import EmbeddingPanel
from tests.test_roi_manager import ImmediateDoc, make_layer, saved, wait_for


class GroupSelectionTests(unittest.TestCase):
    def test_group_visibility_respects_unsequenced_checkbox(self):
        for kind in ('bins', 'spots', 'method_a'):
            with self.subTest(kind=kind):
                layer = make_layer(kind)
                groups = layer.groups_manager
                checkbox = layer.notsequenced_checkbox
                checkbox.active = []
                name = groups.add_mask_group('A', layer.sample.adata.obs_names == 'b0_0',
                                             layer.roi_manager.doc)

                def assert_visible(show_unsequenced, show_group=True):
                    names = list(layer.manager.loaded_obs_names)
                    visible = list(layer.view)
                    self.assertEqual(visible[names.index('b0_0')], show_group)
                    self.assertFalse(visible[names.index('b3_3')])
                    missing = [i for i, value in enumerate(layer.manager.get_nan_vector()) if value]
                    self.assertEqual(len(missing), 1)
                    self.assertEqual(visible[missing[0]], show_unsequenced)

                assert_visible(False)
                checkbox.active = [0]
                assert_visible(True)
                groups._apply_visibility_mask()
                assert_visible(True)
                if layer.lazy_loading:
                    layer.manager.clear_sources()
                    layer.manager.fetch_data()
                    assert_visible(True)
                checkbox.active = []
                groups._apply_visibility_mask()
                assert_visible(False)
                if layer.lazy_loading:
                    layer.manager.clear_sources()
                    layer.manager.fetch_data()
                    assert_visible(False)
                groups._active_groups.discard(name)
                groups._apply_visibility_mask()
                checkbox.active = [0]
                groups._apply_visibility_mask()
                assert_visible(True, show_group=False)

    def test_clear_selection_without_embedding(self):
        layer = make_layer()
        widgets = layer.get_widgets(ImmediateDoc())
        roi = saved(layer.roi_manager)
        layer.roi_manager.select(roi)
        self.assertTrue(layer.get_selected_obs_names())
        clear = next(b for b in widgets.select({'type': Button}) if b.label == 'Clear selection')
        clear._trigger_event(ButtonClick(clear))
        self.assertEqual(layer.get_selected_obs_names(), set())
        self.assertEqual(layer.manager.source.selected.indices, [])

    def test_unselect_and_remove_preserve_unrelated_selection(self):
        for action in ('unselect', 'remove_roi', 'delete_group'):
            with self.subTest(action=action):
                layer = make_layer()
                manager = layer.roi_manager
                roi = saved(manager)
                name = manager.create_group(roi)
                members = set(layer.sample.adata.obs_names[roi.mask])
                layer.set_selected_obs_names(members | {'b3_3'})
                if action == 'delete_group':
                    layer.groups_manager.delete_group(name, manager.doc)
                else:
                    getattr(manager, action)(roi)
                self.assertEqual(layer.get_selected_obs_names(), {'b3_3'})
                displayed = layer.manager.loaded_obs_names
                self.assertEqual({displayed[i] for i in layer.manager.source.selected.indices}, {'b3_3'})

    def test_roi_group_creation_and_removal_log_progress(self):
        layer = make_layer()
        manager = layer.roi_manager
        roi = saved(manager)
        groups = layer.groups_manager
        entered, release = threading.Event(), threading.Event()
        original = groups._prepare_groups

        def blocked(*args):
            entered.set()
            release.wait(5)
            return original(*args)

        with patch.object(groups, '_prepare_groups', side_effect=blocked):
            manager._group_click(roi)
            self.assertTrue(entered.wait(5))
            self.assertTrue(manager.logger.running)
            self.assertIn('Creating group', manager.logger.div.text)
            release.set()
            wait_for(lambda: not groups._operation_running)
        self.assertIn('created', manager.logger.div.text)
        self.assertTrue(any('members' in message for _, message in manager.logger._steps))
        manager.select(roi)
        manager._delete_click(roi)
        wait_for(lambda: not groups._operation_running)
        self.assertIn('deleted', manager.logger.div.text)
        self.assertEqual(layer.get_selected_obs_names(), set())
        self.assertNotIn('ROI: Region', groups._groups)

    def test_worker_defers_updates_to_document_thread(self):
        layer = make_layer()
        groups = layer.groups_manager
        doc = Document()
        groups.get_widgets(doc)
        owner = threading.get_ident()
        workers, commits = [], []

        def work():
            workers.append(threading.get_ident())
            return groups._prepare_groups({'A': [('g0', '>', '15', 'AND')]})

        def done(prepared):
            commits.append(threading.get_ident())
            groups._groups['A'] = [('g0', '>', '15', 'AND')]
            groups._apply_groups_to_adata(prepared)

        groups.run_group_action(doc, 'Creating A', work, done, 'Created A')
        wait_for(lambda: any(isinstance(cb, NextTickCallback) for cb in doc.session_callbacks))
        self.assertNotEqual(workers, [owner])
        self.assertNotIn('grp_A', layer.sample.adata.obs)
        completion = next(cb for cb in doc.session_callbacks if isinstance(cb, NextTickCallback))
        completion.callback()
        self.assertEqual(commits, [owner])
        self.assertIn('grp_A', layer.sample.adata.obs)
        self.assertFalse(groups._operation_running)

    def test_deletion_clears_embedding_selection_too(self):
        layer = make_layer()
        roi = saved(layer.roi_manager)
        group = layer.roi_manager.create_group(roi)
        panel = EmbeddingPanel(ImmediateDoc())
        panel.set_components([layer])
        source = panel._component_views[layer.key]['source']
        names = layer.sample.adata.obs_names.tolist()
        source.data = dict(obs_name=names, x=list(range(16)), y=list(range(16)), color=['red'] * 16)
        layer.roi_manager.select(roi)
        self.assertTrue(source.selected.indices)
        layer.groups_manager.delete_group(group, layer.roi_manager.doc)
        self.assertEqual(source.selected.indices, [])
        self.assertEqual(layer.manager.source.selected.indices, [])

    def test_condition_creation_failure_keeps_conditions_and_no_group(self):
        layer = make_layer()
        groups = layer.groups_manager
        groups._conditions = [('g0', '>', 'invalid', 'AND')]
        widgets = layer.get_widgets(ImmediateDoc())
        create = next(b for b in widgets.select({'type': Button}) if b.label == '✓ Create group')
        create._trigger_event(ButtonClick(create))
        wait_for(lambda: not groups._operation_running)
        self.assertIn('Error', groups.logger.div.text)
        self.assertNotIn('Group A', groups._groups)
        self.assertEqual(groups._conditions, [('g0', '>', 'invalid', 'AND')])

    def test_condition_creation_and_deletion_use_progress_logger(self):
        layer = make_layer()
        groups = layer.groups_manager
        groups._conditions = [('g0', '>', '15', 'AND')]
        widgets = layer.get_widgets(ImmediateDoc())
        buttons = list(widgets.select({'type': Button}))
        create = next(b for b in buttons if b.label == '✓ Create group')
        delete = next(b for b in buttons if b.label == '🗑 Delete group')
        create._trigger_event(ButtonClick(create))
        wait_for(lambda: not groups._operation_running)
        self.assertIn('created', groups.logger.div.text)
        self.assertEqual(int(layer.sample.adata.obs['grp_Group A'].sum()), 8)
        self.assertTrue(any('Evaluating condition' in message for _, message in groups.logger._steps))
        layer.set_selected_obs_names(['b3_3'])
        delete._trigger_event(ButtonClick(delete))
        wait_for(lambda: not groups._operation_running)
        self.assertIn('deleted', groups.logger.div.text)
        self.assertEqual(layer.get_selected_obs_names(), set())


if __name__ == '__main__':
    unittest.main()
