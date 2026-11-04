"""Registration and integration checks for the aligned bin coordinate system."""
import json
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import shapely
from PIL import Image
from bokeh.document import Document
from bokeh.models import Circle, Rect, ColumnDataSource
from bokeh.plotting import figure
from shapely.geometry import Point, box
from shapely.affinity import rotate

from hd_viewer.startup.data.image.image import ImageProcessor
from hd_viewer.startup.data.image.transform import SpatialTransform
from hd_viewer.startup.operations.bins2spots import Bins2Spots
from hd_viewer.startup.sample import Bins, Cells, Spots
from hd_viewer.board.components.layer import Circles, Squares
from hd_viewer.board.utils.measure import MeasureTool
from hd_viewer.board.layout.thumbnail import ViewportReadout as ROIReadout
from hd_viewer.board.operations.layer.de import DEAnalysis
from hd_viewer.board.layout.layer.roi import ROIManager
from hd_viewer.board.components.layer.utils.selection import aggregation_names, native_selection, selection_indices


def fixture(angle=23., reflected=False, skew=.4, image_scale=1., n=6, residual=False):
    theta = np.deg2rad(angle)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    basis = np.array([[4., 0.], [skew, -3.5 if reflected else 3.5]]) @ rotation.T
    rows, cols = np.indices((n, n))
    xy = np.column_stack((cols.ravel(), rows.ravel())) @ basis + [65, 70]
    if residual:
        xy[7] += [.07, -.03]
    names = [f'b{i}' for i in range(len(xy))]
    positions = pd.DataFrame(dict(array_col=cols.ravel(), array_row=rows.ravel(),
         pxl_col_in_fullres=xy[:, 0], pxl_row_in_fullres=xy[:, 1], in_tissue=1, index=names), index=names)
    transform = SpatialTransform.fit(positions, (140, 180, 3), image_scale)
    counts = np.column_stack((np.arange(len(xy)) + 1, np.ones(len(xy))))
    data = ad.AnnData(counts, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['gene', 'one']))
    data.obs['fraction'] = np.arange(len(xy)) / len(xy)
    data.obs['kind'] = pd.Categorical(np.where(cols.ravel() % 2, 'a', 'b'))
    data.uns['mean_obs_columns'] = ['fraction']
    b2s = Bins2Spots(data, positions, pixels_per_bin=4, bins_size=2, spatial_transform=transform)
    sample = Bins(b2s.bin_squares, data, pixels_per_bin=4, bin_size=2, spatial_transform=transform)
    return positions, transform, b2s, sample


class SpatialTransformTests(unittest.TestCase):
    def test_round_trip_bounds_orientation_and_rectangle_registration(self):
        for angle in (0., 2., 179.):
            for reflected in (False, True):
                for scale in (1., .4, 2.):
                    with self.subTest(angle=angle, reflected=reflected, scale=scale):
                        pos, transform, b2s, sample = fixture(angle, reflected, image_scale=scale, residual=True)
                        xy = pos[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy()
                        np.testing.assert_allclose(transform.points(transform.points(xy), inverse=True), xy, atol=1e-10)
                        np.testing.assert_allclose(transform.image_points(xy * scale), transform.points(xy))
                        corners = np.array([[0, 0], [180, 0], [180, 140], [0, 140]]) / scale
                        mapped = transform.points(corners)
                        self.assertGreaterEqual(mapped.min(), -1e-9)
                        self.assertLessEqual(mapped.max(), transform.output_size + 1e-9)
                        # Image coordinates receive no rotation, shear or anisotropic scale.
                        delta = transform.points([[0,0],[1,0],[0,1]])
                        np.testing.assert_allclose(delta[1]-delta[0], [scale,0], atol=1e-10)
                        np.testing.assert_allclose(delta[2]-delta[0], [0,-scale], atol=1e-10)
                        native = sample.fetch_data(colnames=['gene'])
                        for i, row in native.iterrows():
                            rectangle = box(row.x-row.width/2, row.y-row.height/2,
                                            row.x+row.width/2, row.y+row.height/2)
                            rectangle = rotate(rectangle, row.angle, origin=(row.x,row.y), use_radians=True)
                            self.assertLess(rectangle.symmetric_difference(b2s.bin_squares.loc[i].geometry).area, 1e-9)
                        np.testing.assert_allclose(sample.get_coordinates(sample.gdf), transform.points(xy), atol=1e-10)
                        self.assertEqual(list(sample.adata.obs_names), list(pos.index))

    def test_invalid_grids_and_immutable_fingerprint(self):
        positions, transform, _, _ = fixture()
        with self.assertRaises(FrozenInstanceError):
            transform.pitch = 0
        forward = transform.forward
        forward[0, 0] = 0
        self.assertNotEqual(transform.forward[0, 0], 0)
        same = SpatialTransform.fit(positions, (140, 180, 3))
        self.assertEqual(transform.fingerprint, same.fingerprint)
        self.assertNotEqual(transform.fingerprint, SpatialTransform.fit(positions, (140, 180, 3), .5).fingerprint)
        for field in ('array_col', 'pxl_row_in_fullres'):
            broken = positions.copy()
            broken[field] = np.nan
            with self.assertRaisesRegex(ValueError, 'finite'):
                SpatialTransform.fit(broken, (140, 180, 3))
        with self.assertRaisesRegex(ValueError, 'singular'):
            SpatialTransform.fit(positions.iloc[:6], (140, 180, 3))
        with self.assertRaises(ValueError):
            SpatialTransform.fit(positions, (140, 180, 3), 0)

    def test_image_landmarks_and_transparent_padding(self):
        for scale in (1., .5):
            positions, transform, _, _ = fixture(image_scale=scale)
            image = np.zeros((140, 180, 3), np.uint8)
            image[:] = [30, 40, 50]
            xy = positions.iloc[14][['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy(dtype=float)
            x, y = np.rint(xy * scale).astype(int)
            image[y-2:y+3, x-2:x+3] = [240, 5, 20]
            aligned = transform.warp_image(image)
            left = int(transform.forward[0,2])
            top = int(transform.output_size-transform.forward[1,2])
            np.testing.assert_array_equal(aligned[top:top+140,left:left+180,:3], image)
            self.assertTrue((aligned[top:top+140,left:left+180,3]==255).all())
            px, py = transform.points(xy)
            row, col = int(transform.output_size - py), int(px)
            region = aligned[max(0,row-2):row+3, max(0,col-2):col+3]
            self.assertGreater(region[..., 0].max(), 200)
            self.assertEqual(aligned.shape, (transform.output_size, transform.output_size, 4))
            self.assertTrue((aligned[..., 3] == 0).any())

    def test_tiles_and_single_image_share_orientation_and_atomic_caches(self):
        _, transform, _, _ = fixture()
        image = np.random.default_rng(5).integers(0, 256, (140, 180, 3), dtype=np.uint8)
        aligned = transform.warp_image(image)
        with tempfile.TemporaryDirectory() as root:
            processor = ImageProcessor(image, Path(root)/'aligned', tile_size=32, spatial_transform=transform)
            np.testing.assert_array_equal(processor.image, aligned)
            processor.create_tiles(single_image=True)
            with Image.open(processor.tiles_path/'z0_x0_y0.png') as im:
                np.testing.assert_array_equal(im, aligned)
            marker = (processor.tiles_path/'manifest.json').read_text()
            with patch.object(processor, '_write_tiles', side_effect=RuntimeError('interrupted')):
                with self.assertRaises(RuntimeError):
                    processor.create_tiles()
            self.assertEqual((processor.tiles_path/'manifest.json').read_text(), marker)
            self.assertFalse(list(Path(root).glob('.tiles-*')))
            processor.create_tiles()
            for zoom in range(processor.max_zoom + 1):
                count = 2**zoom
                stitched = Image.new('RGBA', (count*32, count*32))
                for x in range(count):
                    for y in range(count):
                        with Image.open(processor.tiles_path/f'z{zoom}_x{x}_y{y}.png') as tile:
                            stitched.paste(tile, (x*32,y*32))
                expected = Image.fromarray(aligned).resize(stitched.size, Image.Resampling.LANCZOS)
                np.testing.assert_array_equal(stitched, expected)
            with patch.object(processor, '_write_tiles', side_effect=AssertionError('cache should be reused')):
                processor.create_tiles(overwrite=False)


class AlignedIntegrationTests(unittest.TestCase):
    def test_array_blocks_are_shared_and_complete_outside_viewport(self):
        _, transform, b2s, sample = fixture(residual=True)
        all_data = sample.fetch_data(colnames=['gene', 'fraction', 'kind'], resolution_level='agg:2')
        mapping = aggregation_names(sample, 2)
        self.assertEqual(mapping.loc['b0'], 'summary:0:0:2')
        self.assertEqual(mapping.loc['b6'], 'summary:0:0:2')
        expected = pd.Series(sample.adata.X[:, 0], index=sample.adata.obs_names).groupby(mapping).sum()
        np.testing.assert_allclose(all_data.gene.reindex(expected.index), expected)
        mean = sample.adata.obs.fraction.groupby(mapping).mean()
        np.testing.assert_allclose(all_data.fraction.reindex(mean.index), mean)
        aggregate = sample._get_aggregated_sample(2).adata
        np.testing.assert_allclose(aggregate.to_df().gene.reindex(expected.index), expected)
        row = all_data.loc['summary:0:0:2']
        # A tiny corner intersects one aggregate but none of its far native bins.
        cropped = sample.fetch_data(bbox=(row.x-row.width/2+.01, row.x-row.width/2+.02,
                                         row.y-row.height/2+.01, row.y-row.height/2+.02),
                                    colnames=['gene'], resolution_level='agg:2')
        self.assertEqual(cropped.loc[row.name, 'gene'], row.gene)
        selected = native_selection(sample, {row.name})
        self.assertEqual(selected, set(mapping.index[mapping==row.name]))
        self.assertEqual(selection_indices(sample, {'b0'}, list(all_data.obs_name)), [list(all_data.obs_name).index(row.name)])
        masks = DEAnalysis(sample).masks_at_level({'group': np.asarray(mapping==row.name)}, 2, aggregate.obs_names)
        self.assertEqual(masks['group'].sum(), 1)
        self.assertTrue(masks['group'][aggregate.obs_names.get_loc(row.name)])
        # Missing native bins and translated image bounds do not redefine blocks.
        sample.gdf = sample.gdf.drop(index='b1')
        sample._grid_cache.clear()
        self.assertEqual(sample.aggregation_names(2).loc['b6'], row.name)
        partial = sample.fetch_data(resolution_level='agg:4', colnames=['gene'])
        self.assertEqual(len(partial), 4)
        self.assertAlmostEqual(partial.gene.sum(), sample.adata.X[:, 0].sum()-sample.adata.X[1, 0])

    def test_startup_builds_aligned_image_before_native_layers(self):
        from hd_viewer.startup import StartupController
        from tests.test_deconvolution import tenx_reference_bytes
        from tests.test_roi_manager import ImmediateDoc
        positions, _, _, _ = fixture()
        positions = positions.rename(index={'b0':'a', 'b1':'b', 'b6':'c'})
        positions['index'] = positions.index
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            image = root/'image.png'
            Image.fromarray(np.full((140,180,3), 180, dtype=np.uint8)).save(image)
            counts = root/'counts.h5'
            counts.write_bytes(tenx_reference_bytes())
            pos_path = root/'positions.parquet'
            positions.drop(columns='index').rename_axis('barcode').reset_index().to_parquet(pos_path)
            (root/'scalefactors_json.json').write_text(json.dumps({'spot_diameter_fullres':4}))
            config = dict(image_path=str(image), expression_path=str(counts),
                          positions_path=str(pos_path), bins_size=2)
            with patch.dict('hd_viewer.startup.AVAILABLE_DATASETS', {'Human CRC':config}), \
                 patch('hd_viewer.startup.data.image.image.ensure_tile_server'):
                startup = StartupController(ImmediateDoc())
                startup.data_source.active = 1
                startup.example_select.value = 'Human CRC'
                startup.pseudospots_checkbox.active = False
                startup.tasks = {'tiles':False, 'bins':False}
                startup.preparing = True
                with patch.object(startup, 'generate_segmentations'):
                    startup._prepare_data()
                self.assertTrue(all(startup.tasks.values()))
                self.assertEqual(len(startup.app.geometries), 1)
                layer = startup.app.geometries[0]
                self.assertIs(layer.sample.spatial_transform, startup.spatial_transform)
                self.assertEqual(layer.sample.adata.n_obs, 3)
                self.assertTrue((startup.image_processor.tiles_path/'manifest.json').exists())
                self.assertIn('padded-v2-', str(startup.image_processor.tiles_path))
                layer.get_widgets(startup.doc)
                layer.refresh_glyph()
                self.assertIsInstance(layer.glyph_renderer.glyph, Rect)
                self.assertNotIn('xs', layer.manager.source.data)
                self.assertNotIn('ys', layer.manager.source.data)
                startup.doc.add_root(startup.app.figure_panel.p)
                startup.doc.to_json()
                # Reopen the same example: neither tile reuse nor bin loading
                # should materialize the original or padded image.
                with patch('hd_viewer.startup.data.VisiumHD.VisiumHD.read_img', side_effect=AssertionError('Cached example decoded pixels')), \
                     patch.object(SpatialTransform, 'warp_image', side_effect=AssertionError('Cached example padded pixels')):
                    reopened = StartupController(ImmediateDoc())
                    reopened.data_source.active = 1
                    reopened.example_select.value = 'Human CRC'
                    reopened.pseudospots_checkbox.active = False
                    reopened.tasks = {'tiles':False, 'bins':False}
                    reopened.preparing = True
                    with patch.object(reopened, 'generate_segmentations'):
                        reopened._prepare_data()
                    self.assertTrue(all(reopened.tasks.values()))
                    self.assertIsNone(reopened.dataset.img)
                    self.assertIsNone(reopened.image_processor._image)
                    self.assertEqual(len(reopened.app.geometries), 1)

    def test_metric_radii_readout_ruler_and_polygon_layers(self):
        _, transform, b2s, sample = fixture(skew=1.4)
        plot = figure(x_range=(20, 60), y_range=(30, 80))
        component = Squares(plot, sample, lazy_loading=True, alias='Bin')
        manager = component.roi_manager
        readout = ROIReadout(plot)
        readout.set_transform(transform, .5)
        manager.coordinate_scale = readout.scale
        factors, unit = manager.radius_scale()
        self.assertEqual(unit, 'µm')
        seed = transform.points([70., 75.])
        shape = manager.seed_radius_geometry(*seed, 5., *factors)
        original = transform.geometry(shape, inverse=True)
        np.testing.assert_allclose(original.bounds, [60,65,80,85], atol=1e-9)
        sample.adata.obs['grp_seeds'] = False
        sample.adata.obs.loc[['b0','b35'], 'grp_seeds'] = True
        combined, count = manager.compute_group_radius('seeds', 5, factors)
        self.assertEqual(count, 2)
        centers = sample.get_coordinates(sample.gdf.loc[['b0','b35']]).to_numpy()
        expected = shapely.union_all([manager.seed_radius_geometry(*xy,5,*factors) for xy in centers])
        self.assertLess(expected.symmetric_difference(combined).area, 1e-8)
        corners = transform.points([[20,30],[60,30],[60,80],[20,80]], inverse=True)
        self.assertEqual(readout.bbox(), tuple(round(v) for v in (*corners.min(axis=0), *corners.max(axis=0))))
        measure = MeasureTool(plot)
        measure.set_transform(transform, .5)
        s = measure.scale.data
        metric = np.array([[s['a'][0],s['b'][0]],[s['c'][0],s['d'][0]]])
        delta = transform.points([80.,90.])-transform.points([70.,75.])
        self.assertAlmostEqual(np.linalg.norm(metric @ delta), np.hypot(10,15)*.5)
        # Segmentation coordinates are in image pixels, not full-res pixels.
        cell = Point(70,75).buffer(2)
        mapped = transform.image_geometry(cell)
        self.assertTrue(mapped.equals_exact(transform.geometry(cell), 1e-9))
        _, _, spot_b2s, _ = fixture(n=60)
        self.assertGreater(len(spot_b2s.spot_polygons), 0)
        spots = Spots(spot_b2s.spot_polygons, ad.AnnData(np.ones((len(spot_b2s.spot_polygons),1)),
                       obs=pd.DataFrame(index=spot_b2s.spot_polygons.index.astype(str)), var=pd.DataFrame(index=['g'])))
        spots.spatial_transform = spot_b2s.spatial_transform
        source = spots.fetch_data()
        self.assertTrue({'x', 'y', 'radius'}.issubset(source.columns))
        self.assertNotIn('xs', source)
        self.assertNotIn('ys', source)
        self.assertTrue(all(np.isscalar(x) for x in source.x))
        self.assertTrue(all(np.isscalar(y) for y in source.y))
        expected_radius = spot_b2s.spatial_transform.pitch * 55. / (2 * spot_b2s.bins_size)
        np.testing.assert_allclose(source.radius, expected_radius)
        circles = Circles('radius', plot, spots, alias='Spot')
        circles.manager.fetch_data()
        circles.add_glyph()
        self.assertIsInstance(circles.glyph_renderer.glyph, Circle)
        self.assertEqual(circles.glyph_renderer.glyph.radius, 'radius')
        component.manager.clear_sources()
        component.add_glyph()
        self.assertIsInstance(component.glyph_renderer.glyph, Rect)
        doc = Document()
        doc.add_root(plot)
        doc.validate()
        doc.to_json()


if __name__ == '__main__':
    unittest.main()
