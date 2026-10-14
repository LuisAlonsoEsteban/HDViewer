import math
import unittest
from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from shapely.geometry import Polygon

from bokeh_server.components.bin_info import BinInfoOverlay
from bokeh_server.loader import LazyLoader
from bokeh_server.manager import DynamicSourceManager
from src.sample import Bins


class BinInfoOverlayTests(unittest.TestCase):
    def setUp(self):
        # A native bin transformed by rotation followed by unequal axis scaling.
        angle = math.radians(30)
        rotation = np.array([[math.cos(angle), -math.sin(angle)],
                             [math.sin(angle), math.cos(angle)]])
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) @ rotation.T
        corners = corners * [1, 2] + [10, 10]
        self.corners = corners
        gdf = gpd.GeoDataFrame({'index': ['bin']}, geometry=[Polygon(corners)],
                               index=pd.Index(['bin'], name='barcode'))
        adata = ad.AnnData(np.ones((1, 1)), obs=pd.DataFrame(index=['bin']),
                          var=pd.DataFrame(index=['gene']))
        self.sample = Bins(gdf, adata, pixels_per_bin=2, tile_scale_yfactor=2,
                           native_angle=89, bin_size_um=16)
        self.manager = DynamicSourceManager(self.sample)
        self.plot = figure(x_range=(0, 20), y_range=(0, 20))
        self.overlay = BinInfoOverlay(SimpleNamespace(
            sample=self.sample, manager=self.manager, figure=self.plot))

    def test_native_dimensions_and_rotation_use_polygon_edges(self):
        width, height, angle = self.overlay.metadata()['geometry']
        edge = self.corners[1] - self.corners[0]
        np.testing.assert_allclose(width, np.linalg.norm(edge))
        np.testing.assert_allclose(height, np.linalg.norm(self.corners[2] - self.corners[1]))
        self.assertAlmostEqual(angle, math.degrees(math.atan2(edge[1], edge[0])))
        self.assertIn('1 × 1 · native', self.overlay.header.text)
        self.assertIn('Native bin 16 µm', self.overlay.label.text)
        self.assertIn('plot units', self.overlay.label.text)

    def test_zoom_updates_overlay_via_source_callbacks(self):
        loader = LazyLoader(self.manager, self.plot)
        self.manager.clear_sources()
        self.plot.x_range.end = 10000
        self.plot.y_range.end = 10000
        loader.force_reload()
        self.assertEqual(self.overlay.metadata()['factor'], 12)
        self.assertIn('12 × 12', self.overlay.header.text)
        self.assertIn('Rotation   0.00°', self.overlay.label.text)
        self.assertEqual(self.overlay.metadata()['geometry'], (24, 48, 0))
        self.plot.x_range.end = 20
        self.plot.y_range.end = 20
        loader.force_reload()
        self.assertIn('1 × 1 · native', self.overlay.header.text)

    def test_cluster_locked_factor_overrides_zoom(self):
        self.sample.adata.obs['leiden_4'] = ['0']
        self.manager.get_colnames = lambda: ['leiden_4']
        self.manager.set_fetch_context(resolution_level='agg:12')
        self.manager.clear_sources()
        self.assertEqual(self.overlay.metadata()['factor'], 4)
        self.assertIn('4 × 4', self.overlay.header.text)
        self.assertIn('Cluster-locked', self.overlay.label.text)
        self.manager.get_colnames = lambda: []
        self.manager.clear_sources()
        self.assertIn('12 × 12', self.overlay.header.text)
        self.assertIn('Zoom-driven', self.overlay.label.text)

    def test_screen_anchor_visibility_and_renderer_replacement(self):
        label = self.overlay.panel
        self.assertEqual((label.position.target, label.position.symbol), ('frame', 'top_right'))
        self.assertEqual(label.anchor, 'top_right')
        source = ColumnDataSource({'xs': [], 'ys': []})
        first = self.plot.patches(xs='xs', ys='ys', source=source)
        self.overlay.bind_renderer(first)
        first.visible = False
        self.assertFalse(label.visible)
        second = self.plot.patches(xs='xs', ys='ys', source=source)
        self.overlay.bind_renderer(second)
        first.visible = True
        second.visible = False
        first.visible = False
        first.visible = True
        self.assertFalse(label.visible)
        second.visible = True
        self.assertTrue(label.visible)
        doc = Document()
        doc.add_root(self.plot)
        doc.to_json()
        self.assertEqual(len(list(self.plot.select({'name': 'bin_resolution_info'}))), 1)

    def test_fixed_level_survives_zoom_then_dynamic_uses_current_view(self):
        self.overlay.mode_select.value = '4'
        self.assertEqual(self.manager.fixed_aggregation_factor, 4)
        self.assertEqual(self.overlay.metadata()['factor'], 4)
        self.assertIn('Resolution Fixed', self.overlay.label.text)
        self.assertTrue(all(str(name).endswith(':4') for name in self.manager.source.data['obs_name']))
        loader = LazyLoader(self.manager, self.plot)
        self.plot.x_range.end = 10000
        self.plot.y_range.end = 10000
        loader.force_reload()
        self.assertEqual(self.overlay.metadata()['factor'], 4)
        self.overlay.mode_select.value = 'dynamic'
        self.assertIsNone(self.manager.fixed_aggregation_factor)
        self.assertEqual(self.overlay.metadata()['factor'], 12)
        self.assertIn('Zoom-driven', self.overlay.label.text)

    def test_raw_mode_overrides_zoom_and_reloads_native_bins(self):
        self.plot.x_range.end = 10000
        self.plot.y_range.end = 10000
        self.overlay.mode_select.value = '4'
        self.overlay.mode_select.value = '1'
        self.assertEqual(self.manager.fixed_aggregation_factor, 1)
        self.assertIsNone(self.manager.effective_resolution_level)
        self.assertEqual(self.manager.source.data['obs_name'], ['bin'])
        np.testing.assert_allclose(self.overlay.metadata()['geometry'][0],
                                   np.linalg.norm(self.corners[1] - self.corners[0]))
        self.assertIn('1 × 1 · native', self.overlay.header.text)

    def test_cluster_lock_restores_saved_mode(self):
        self.overlay.mode_select.value = '2'
        self.sample.adata.obs['leiden_4'] = ['0']
        self.manager.get_colnames = lambda: ['leiden_4']
        self.manager.clear_sources()
        self.assertTrue(self.overlay.mode_select.disabled)
        self.assertTrue(self.overlay.lock_notice.visible)
        self.assertEqual(self.overlay.mode_select.value, '4')
        self.assertEqual(self.manager.fixed_aggregation_factor, 2)
        # Even programmatic changes while disabled cannot silently change mode.
        self.overlay.mode_select.value = '1'
        self.assertEqual(self.overlay.mode_select.value, '4')
        self.assertEqual(self.manager.fixed_aggregation_factor, 2)
        self.manager.get_colnames = lambda: []
        self.manager.clear_sources()
        self.assertFalse(self.overlay.mode_select.disabled)
        self.assertFalse(self.overlay.lock_notice.visible)
        self.assertEqual(self.overlay.mode_select.value, '2')
        self.assertEqual(self.overlay.metadata()['factor'], 2)

    def test_supported_fixed_levels_and_invalid_factor(self):
        self.assertEqual([value for value, _ in self.overlay.mode_select.options],
                         ['dynamic', '1', '2', '4', '8', '12'])
        for factor in self.sample.AGGREGATION_FACTORS:
            self.overlay.mode_select.value = str(factor)
            self.assertEqual(self.overlay.metadata()['factor'], factor)
        with self.assertRaises(ValueError):
            self.manager.set_aggregation_factor(3)
        self.assertEqual(self.manager.fixed_aggregation_factor, 12)

    def test_empty_geometry(self):
        self.sample.gdf = self.sample.gdf.iloc[:0]
        overlay = BinInfoOverlay(SimpleNamespace(
            sample=self.sample, manager=self.manager, figure=figure()))
        self.assertIn('No bin geometry available', overlay.label.text)

    def test_minimize_keeps_the_current_resolution_visible(self):
        self.assertFalse(self.overlay.minimized)
        self.assertTrue(self.overlay.body.visible)

        self.overlay.toggle_button._trigger_event(ButtonClick(self.overlay.toggle_button))
        self.assertTrue(self.overlay.minimized)
        self.assertFalse(self.overlay.body.visible)
        self.assertEqual(self.overlay.toggle_button.label, '+')
        # The header still reports the active resolution, and keeps updating
        self.assertIn('1 × 1 · native', self.overlay.header.text)
        self.overlay.mode_select.value = '4'
        self.assertIn('4 × 4', self.overlay.header.text)

        self.overlay.toggle_button._trigger_event(ButtonClick(self.overlay.toggle_button))
        self.assertFalse(self.overlay.minimized)
        self.assertTrue(self.overlay.body.visible)
        self.assertEqual(self.overlay.toggle_button.label, '−')


if __name__ == '__main__':
    unittest.main()
