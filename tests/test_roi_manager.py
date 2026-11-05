import threading
import time
import unittest
from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick, SelectionGeometry, Tap
from bokeh.layouts import column
from bokeh.models import Button, ColumnDataSource
from bokeh.plotting import figure
from shapely.geometry import box

from hd_viewer.board.components.layer import Circles, Patches, Squares
from hd_viewer.board.layout.layer.roi import ROI, ROIManager
from hd_viewer.board.components.layer.utils.selection import aggregation_names
from hd_viewer.board.layout.embedding import EmbeddingPanel
from hd_viewer.startup.sample import Bins, Cells, Spots


class ImmediateDoc(Document):
    def add_next_tick_callback(self, callback):
        callback()

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass


def make_layer(kind='bins', plot=None):
    names = [f'b{x}_{y}' for y in range(4) for x in range(4)]
    all_names = names + ['not_sequenced']
    polygons = [box(x, y, x + 1, y + 1) for y in range(4) for x in range(4)] + [box(1.2, 1.2, 1.4, 1.4)]
    gdf = gpd.GeoDataFrame({'index': all_names}, geometry=polygons,
                           index=pd.Index(all_names, name='barcode'))
    counts = np.arange(32, dtype=float).reshape(16, 2) + 1
    data = ad.AnnData(counts, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['g0', 'g1']))
    plot = plot or figure(x_range=(-1, 5), y_range=(-1, 5))
    if kind == 'bins':
        sample = Bins(gdf, data, pixels_per_bin=1, bin_size=2)
        component = Squares(plot, sample, lazy_loading=True, alias='Bin')
    elif kind == 'spots':
        component = Circles(.5, plot, Spots(gdf, data), alias='Spot')
    else:
        component = Patches(plot, Cells(gdf, data), lazy_loading=True, alias='Cell', key=f'cells_{kind}')
    component.add_parameter(attr='alpha', colname='g0', colname_hovered=True)
    component.add_parameter(attr='color', static=True, val='#6baed6', not_sequenced_val='#b8b8b8')
    if component.lazy_loading:
        component.manager.bbox = (-1, 5, -1, 5)
    component.manager.fetch_data()
    if component.lazy_loading:
        component.loader = SimpleNamespace(force_reload=component.manager.fetch_data)
    component.get_widgets(ImmediateDoc())
    component._figure_de_layout = column()
    return component


def saved(manager, name='Region', geometry=None):
    geometry = box(.2, .2, 1.8, 1.8) if geometry is None else geometry
    roi = ROI(name, geometry, 'box', '#ff0000', mask=manager.compute_membership(geometry))
    manager.add_roi(roi)
    return roi


def wait_for(predicate, message='Timed out'):
    deadline = time.monotonic() + 15
    while not predicate() and time.monotonic() < deadline:
        time.sleep(.01)
    if not predicate():
        raise AssertionError(message)


class ROIGeometryTests(unittest.TestCase):
    def test_membership_uses_intersection_and_full_native_data(self):
        for kind in ('spots', 'bins', 'method_a'):
            with self.subTest(kind=kind):
                layer = make_layer(kind)
                manager = layer.roi_manager
                self.assertIn('ROIs', [tab.title for tab in layer.get_widgets(ImmediateDoc()).tabs])
                # Nine native shapes touch or intersect the box; an unsequenced shape is excluded.
                mask = manager.compute_membership(box(1, 1, 2, 2))
                self.assertEqual(mask.shape, (16,))
                self.assertEqual(int(mask.sum()), 9)
                if layer.lazy_loading:
                    layer.manager.bbox = (0, 1, 0, 1)
                    layer.manager.fetch_data()
                far = saved(manager, 'Far', box(2.1, 2.1, 3.9, 3.9))
                self.assertEqual(far.count, 4)
                self.assertEqual(manager.compute_membership(box(10, 10, 11, 11)).sum(), 0)
                lasso = manager.lasso_geometry([.2, 1.8, 1.8, .2], [.2, .2, 1.8, 1.8])
                np.testing.assert_array_equal(manager.compute_membership(lasso), saved(manager).mask)

    def test_geometry_validation_and_unequal_radius_scales(self):
        invalid = [lambda: ROIManager.box_geometry(0, 0, 0, 1),
                   lambda: ROIManager.box_geometry(0, 0, np.inf, 1),
                   lambda: ROIManager.lasso_geometry([0, 1], [0, 1]),
                   lambda: ROIManager.lasso_geometry([0, 1, 0, 1], [0, 1, 1, 0]),
                   lambda: ROIManager.lasso_geometry([0, 1, np.nan], [0, 1, 2]),
                   lambda: ROIManager.radius_geometry(0, 0, 0),
                   lambda: ROIManager.radius_geometry(0, 0, -1)]
        for draw in invalid:
            with self.assertRaises(ValueError):
                draw()
        shape = ROIManager.radius_geometry(5, 6, 2, 2, .5)
        np.testing.assert_allclose(shape.bounds, (4, 2, 6, 10))
        np.testing.assert_allclose(ROIManager.box_geometry(2, 3, 0, 1).bounds, (0, 1, 2, 3))

    def test_focus_visibility_and_empty_region(self):
        layer = make_layer()
        manager = layer.roi_manager
        roi = saved(manager)
        old_x, old_y = layer.figure.x_range, layer.figure.y_range
        manager.focus(roi)
        self.assertIs(layer.figure.x_range, old_x)
        self.assertIs(layer.figure.y_range, old_y)
        np.testing.assert_allclose([old_x.start, old_y.start, old_x.end, old_y.end], [.12, .12, 1.88, 1.88])
        manager._toggle_outline(roi, False)
        self.assertEqual(manager.source.data['xs'], [])
        manager._toggle_outline(roi, True)
        self.assertEqual(len(manager.source.data['xs']), 1)
        empty = saved(manager, 'Empty', box(10, 10, 11, 11))
        self.assertEqual(empty.count, 0)
        controls = manager.roi_list.children[-1]
        buttons = list(controls.select({'type': Button}))
        self.assertTrue(next(b for b in buttons if b.label == 'Select').disabled)
        self.assertTrue(next(b for b in buttons if b.label == 'Create group').disabled)
        self.assertFalse(next(b for b in buttons if b.label == 'Focus').disabled)
        manager.select(empty)
        self.assertIsNone(manager.create_group(empty))
        manager.remove_roi(empty.id)
        self.assertNotIn(empty, manager.get_rois())


class ROIDrawingTests(unittest.TestCase):
    def test_deleted_pending_roi_and_failed_membership_do_not_leave_active_controls(self):
        layer = make_layer()
        manager = layer.roi_manager
        real = manager.compute_membership
        running, release = threading.Event(), threading.Event()
        def blocked(geometry):
            running.set()
            release.wait(10)
            return real(geometry)
        manager.compute_membership = blocked
        manager.set_draft(box(.2, .2, 1.8, 1.8), 'box')
        manager.save_roi()
        self.assertTrue(running.wait(5))
        manager.remove_roi(manager.get_rois()[0])
        release.set()
        wait_for(lambda: not manager.logger.running)
        self.assertEqual(manager.get_rois(), [])
        self.assertEqual(manager.source.data['xs'], [])
        self.assertIn('deleted', manager.logger.div.text)

        def failing(geometry):
            raise ValueError('Membership unavailable')
        manager.compute_membership = failing
        manager.set_draft(box(.2, .2, 1.8, 1.8), 'box')
        manager.save_roi()
        wait_for(lambda: manager.get_rois()[0].error is not None)
        self.assertIn('Membership unavailable', manager.logger.div.text)
        self.assertTrue(next(b for b in manager.roi_list.select({'type': Button}) if b.label == 'Select').disabled)

    def test_exclusive_tools_final_geometry_and_selection_are_independent(self):
        plot = figure(x_range=(0, 5), y_range=(0, 5))
        first = make_layer('bins', plot)
        second = make_layer('method_a', plot)
        one, two = first.roi_manager, second.roi_manager
        initial_tools = len(plot.tools)
        first.get_widgets(ImmediateDoc())
        self.assertEqual(len(plot.tools), initial_tools)
        self.assertIs(one.drawing, two.drawing)
        self.assertEqual(len(plot._event_callbacks['selectiongeometry']), 1)
        first.set_selected_obs_names(['b0_0'])
        previous = plot.toolbar.active_drag, plot.toolbar.active_tap
        one.begin()
        self.assertEqual(one.drawing.box_tool.renderers, [])
        plot._trigger_event(SelectionGeometry(plot, geometry={'type': 'rect', 'x0': .2, 'y0': .2, 'x1': 1.8, 'y1': 1.8}, final=False))
        self.assertIsNone(one.draft)
        two.begin()
        self.assertIs(one.drawing.owner, two)
        plot._trigger_event(SelectionGeometry(plot, geometry={'type': 'rect', 'x0': .2, 'y0': .2, 'x1': 1.8, 'y1': 1.8}, final=True))
        self.assertIsNone(one.draft)
        self.assertIsNotNone(two.draft)
        self.assertEqual((plot.toolbar.active_drag, plot.toolbar.active_tap), previous)
        self.assertEqual(first.get_selected_obs_names(), {'b0_0'})
        self.assertEqual(second.get_selected_obs_names(), set())
        two.cancel()
        two.method.value = 'lasso'
        two.begin()
        plot._trigger_event(SelectionGeometry(plot, geometry={'type': 'poly', 'x': [.2, 1.8, 1.8, .2], 'y': [.2, .2, 1.8, 1.8]}))
        self.assertEqual(two.draft_kind, 'lasso')
        self.assertEqual(two.get_rois(), [])

    def test_seed_radius_uses_calibration_or_full_resolution_pixels(self):
        layer = make_layer()
        manager = layer.roi_manager
        manager.coordinate_scale = ColumnDataSource(dict(px_x=[2.], px_y=[1.], um_per_px=[.5]))
        manager.method.value = 'radius'
        manager.radius.value = 2.
        manager.begin()
        self.assertIn('µm', manager.radius.title)
        layer.figure._trigger_event(Tap(layer.figure, x=2, y=2, sx=0, sy=0))
        np.testing.assert_allclose(manager.draft.bounds, (0, -2, 4, 6))
        manager.cancel()
        manager.coordinate_scale.data = dict(px_x=[2.], px_y=[4.], um_per_px=[None])
        manager.begin()
        self.assertIn('full-res px', manager.radius.title)
        layer.figure._trigger_event(Tap(layer.figure, x=2, y=2, sx=0, sy=0))
        np.testing.assert_allclose(manager.draft.bounds, (1, 1.5, 3, 2.5))
        manager.cancel()
        manager.radius.value = -1
        manager.begin()
        self.assertIsNone(manager.drawing.owner)
        self.assertIn('positive radius', manager.logger.div.text)

    def test_save_computes_membership_in_worker_and_rebuilds_without_duplicate_models(self):
        layer = make_layer()
        manager = layer.roi_manager
        real = manager.compute_membership
        computing, release = threading.Event(), threading.Event()
        def slow(geometry):
            computing.set()
            release.wait(10)
            return real(geometry)
        manager.compute_membership = slow
        manager.set_draft(manager.box_geometry(.2, .2, 1.8, 1.8), 'box')
        manager.name_input.value = 'Saved'
        manager.save_button._trigger_event(ButtonClick(manager.save_button))
        self.assertTrue(computing.wait(5))
        roi = manager.get_rois()[0]
        self.assertIsNone(roi.mask)
        self.assertTrue(next(b for b in manager.roi_list.select({'type': Button}) if b.label == 'Select').disabled)
        release.set()
        wait_for(lambda: roi.mask is not None and not manager.logger.running)
        self.assertEqual(roi.count, 4)
        self.assertIs(manager.get_widgets(ImmediateDoc())[0], manager._layout)
        document = Document()
        document.add_root(column(layer.figure, layer.get_widgets(ImmediateDoc())))
        document.validate()
        document.to_json()
        manager.set_draft(box(.1, .1, .9, .9), 'box')
        manager.name_input.value = 'Saved'
        manager.save_roi()
        self.assertEqual(len(manager.get_rois()), 1)
        self.assertIn('already exists', manager.logger.div.text)
        manager.set_visible(False)
        self.assertFalse(manager.renderer.visible)
        self.assertIsNone(manager.draft)
        manager.dispose()
        self.assertNotIn(manager.renderer, layer.figure.renderers)


class ROIGroupTests(unittest.TestCase):
    def test_de_in_flight_is_discarded_when_linked_roi_is_deleted(self):
        layer = make_layer()
        manager, groups = layer.roi_manager, layer.groups_manager
        roi = saved(manager)
        name = manager.create_group(roi)
        real = layer.de_manager.de_analysis.run
        running, release = threading.Event(), threading.Event()
        def blocked(*args, **kwargs):
            running.set()
            release.wait(10)
            return real(*args, **kwargs)
        layer.de_manager.de_analysis.run = blocked
        layer.de_manager._run_de_thread([name], 0, manager.doc)
        self.assertTrue(running.wait(5))
        manager.remove_roi(roi)
        release.set()
        wait_for(lambda: not layer.de_manager._de_logger.running)
        self.assertEqual(layer.de_manager._de_results, {})
        self.assertEqual(layer._figure_de_layout.children, [])
        self.assertIn('discarded', layer.de_manager._de_status_div.text)

    def test_fixed_groups_survive_conditions_and_delete_with_roi(self):
        layer = make_layer()
        manager, groups = layer.roi_manager, layer.groups_manager
        roi = saved(manager)
        group_name = manager.create_group(roi)
        np.testing.assert_array_equal(layer.sample.adata.obs[f'grp_{group_name}'], roi.mask)
        groups._groups['Condition'] = [('g0', '>', '15', 'AND')]
        groups._active_groups.add('Condition')
        groups._apply_groups_to_adata()
        np.testing.assert_array_equal(layer.sample.adata.obs[f'grp_{group_name}'], roi.mask)
        self.assertIn(group_name, layer.de_manager._de_options('groups'))
        layer.de_manager._de_clusters = [group_name]
        layer.de_manager._de_results_source = 'groups'
        layer._figure_de_layout.children = [column()]
        manager.remove_roi(roi)
        self.assertNotIn(group_name, groups._groups)
        self.assertNotIn(f'grp_{group_name}', layer.sample.adata.obs)
        self.assertEqual(layer._figure_de_layout.children, [])
        self.assertIn('Condition', groups._groups)

    def test_manual_group_delete_recreation_and_unique_names(self):
        layer = make_layer()
        manager, groups = layer.roi_manager, layer.groups_manager
        roi = saved(manager)
        groups.add_mask_group('ROI: Region', np.arange(16) < 4, manager.doc)
        linked = manager.create_group(roi)
        self.assertEqual(linked, 'ROI: Region (2)')
        self.assertEqual(manager.create_group(roi), linked)
        groups.delete_group(linked, manager.doc)
        self.assertIsNone(groups.group_for_roi(roi.id))
        self.assertFalse(next(b for b in manager.roi_list.select({'type': Button}) if b.label == 'Create group').disabled)
        recreated = manager.create_group(roi)
        self.assertEqual(recreated, linked)
        manager.remove_roi(roi)
        self.assertIn('ROI: Region', groups._groups)
        groups.delete_group('ROI: Region', manager.doc)
        self.assertEqual(groups._groups, {})
        self.assertNotIn('custom_groups', layer.sample.adata.obs)
        self.assertTrue(layer.manager.get_param('color').is_static)

    def test_roi_group_runs_de_through_existing_controls(self):
        layer = make_layer()
        manager, groups = layer.roi_manager, layer.groups_manager
        roi = saved(manager)
        name = manager.create_group(roi)
        layer.de_manager._run_de_thread([name], 0, manager.doc)
        # `running` clears before on_done has attached the figures: wait for those.
        wait_for(lambda: layer._figure_de_layout.children, layer.de_manager._de_status_div.text)
        self.assertIn(name, layer.de_manager._de_results, layer.de_manager._de_status_div.text)
        manager.remove_roi(roi)
        self.assertEqual(layer._figure_de_layout.children, [])


class ROISelectionTests(unittest.TestCase):
    def test_native_selection_projects_to_aggregated_map_and_embedding(self):
        layer = make_layer()
        roi = saved(layer.roi_manager, geometry=box(.1, .1, .9, .9))
        panel = EmbeddingPanel(ImmediateDoc())
        panel.set_components([layer])
        view = panel._component_views[layer.key]
        layer.manager.set_fetch_context(bbox=(-1, 5, -1, 5), resolution_level='agg:2')
        layer.manager.fetch_data()
        names = list(layer.manager.loaded_obs_names)
        view['source'].data = dict(obs_name=names, x=list(range(len(names))), y=list(range(len(names))), color=['red'] * len(names))
        layer.roi_manager.select(roi)
        self.assertEqual(layer.get_selected_obs_names(), {'b0_0'})
        expected = names.index('summary:0:0:2')
        self.assertEqual(layer.manager.source.selected.indices, [expected])
        self.assertEqual(view['source'].selected.indices, [expected])
        # Selecting an aggregate expands to its native members; native embeddings see all of them.
        view['source'].selected.indices = [expected]
        panel._on_embedding_selection_change(layer.key, [expected])
        self.assertEqual(layer.get_selected_obs_names(), {'b0_0', 'b1_0', 'b0_1', 'b1_1'})
        native = layer.sample.adata.obs_names.tolist()
        view['source'].data = dict(obs_name=native, x=list(range(16)), y=list(range(16)), color=['red'] * 16)
        panel._update_embedding_selection_from_component(layer.key)
        self.assertEqual(view['source'].selected.indices, [0, 1, 4, 5])
        # Projection did not need expression aggregation.
        self.assertEqual(layer.sample._aggregated_samples, {})
        labels = aggregation_names(layer.sample, 2)
        self.assertEqual(labels.loc['b0_0'], 'summary:0:0:2')
        layer.manager.set_fetch_context(bbox=(-1, 5, -1, 5), resolution_level=None)
        layer.manager.fetch_data()
        self.assertEqual(layer.manager.source.selected.indices, [0, 1, 4, 5])


if __name__ == '__main__':
    unittest.main()
