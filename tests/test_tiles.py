import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
from PIL import Image

from hd_viewer.board.components.image import ImageComponent
from hd_viewer.startup.data.image.image import ImageProcessor, make_tile_app, register_tiles


def tile_dir(root, name, colour):
    path = Path(root) / name
    path.mkdir()
    Image.fromarray(np.full((4, 4, 3), colour, np.uint8)).save(path / "z0_x0_y0.png")
    return path


class TileServerTests(unittest.TestCase):
    def test_each_dataset_is_served_under_its_own_key(self):
        with tempfile.TemporaryDirectory() as root:
            brain = register_tiles(tile_dir(root, "brain", 10))
            ovary = register_tiles(tile_dir(root, "ovary", 200))
            self.assertNotEqual(brain, ovary)
            client = make_tile_app().test_client()
            pixel = lambda key: np.asarray(Image.open(
                __import__("io").BytesIO(client.get(f"/{key}/0/0/0").data)))[0, 0, 0]
            self.assertEqual(pixel(brain), 10)
            self.assertEqual(pixel(ovary), 200)      # not the first dataset's image
            self.assertEqual(client.get("/unknown/0/0/0").status_code, 404)
            self.assertEqual(client.get(f"/{ovary}/3/0/0").status_code, 404)

    def test_image_urls_carry_the_key(self):
        with tempfile.TemporaryDirectory() as root:
            processor = ImageProcessor(np.zeros((8, 8, 3), np.uint8), tile_dir(root, "img", 0))
            processor.deploy_server = lambda port, host=None: register_tiles(processor.tiles_path)   # no real server
            component = ImageComponent(processor)
            component.deploy_server()
            plot = SimpleNamespace(calls=[], image_url=lambda **kw: plot.calls.append(kw))
            component.add_lowres_image(plot)
            self.assertEqual(plot.calls[0]["url"], [f"http://localhost:4000/{component.tile_key}/0/0/0"])


class TileGenerationTests(unittest.TestCase):
    def test_every_zoom_covers_rectangular_and_small_images(self):
        for height, width in ((3, 130), (130, 3), (71, 101), (5, 7)):
            with self.subTest(shape=(height, width)), tempfile.TemporaryDirectory() as root:
                image = np.random.default_rng(0).integers(0, 256, (height, width, 3), dtype=np.uint8)
                processor = ImageProcessor(image, Path(root) / "tiles", tile_size=16)
                progress = []
                processor.create_tiles(progress_callback=lambda done, total: progress.append((done, total)))
                for zoom in range(processor.max_zoom + 1):
                    count = 2 ** zoom
                    size = processor.tile_size * count
                    stitched = Image.new("RGB", (size, size))
                    for x in range(count):
                        for y in range(count):
                            path = processor.tiles_path / processor.filename_pattern.format(zoom=zoom, x=x, y=y)
                            with Image.open(path) as tile:
                                self.assertEqual(tile.size, (16, 16))
                                stitched.paste(tile, (x * 16, y * 16))
                    expected = Image.fromarray(np.flipud(image)).resize((size, size), Image.LANCZOS)
                    np.testing.assert_array_equal(np.asarray(stitched), np.asarray(expected))
                total = sum(4 ** zoom for zoom in range(processor.max_zoom + 1))
                self.assertEqual(progress, [(done, total) for done in range(1, total + 1)])

    def test_cache_hit_never_loads_or_pads_image_pixels(self):
        from tests.test_spatial_alignment import fixture
        from hd_viewer.startup.data.image.transform import SpatialTransform
        _, transform, _, _ = fixture()
        image = np.zeros((140,180,3), np.uint8)
        with tempfile.TemporaryDirectory() as root:
            load = Mock(return_value=image)
            first = ImageProcessor(None, Path(root)/'tiles', spatial_transform=transform, image_loader=load)
            self.assertEqual(first.shape, (transform.output_size, transform.output_size, 4))
            self.assertGreater(first.max_zoom, -1)
            self.assertGreater(first.mem, 0)
            self.assertEqual(first.scale_xfactor, 1)
            load.assert_not_called()
            self.assertTrue(first.create_tiles(overwrite=False))
            load.assert_called_once_with()
            cached = ImageProcessor(None, Path(root)/'tiles', spatial_transform=transform,
                                    image_loader=Mock(side_effect=AssertionError('Cache must not decode pixels')))
            with patch.object(SpatialTransform, 'warp_image', side_effect=AssertionError('Cache must not pad pixels')):
                self.assertFalse(cached.create_tiles(overwrite=False))
                self.assertIsNone(cached._image)
                self.assertFalse(ImageComponent(cached).are_tiles_needed)
            # A missing manifest forces a fresh generation and exactly one decode.
            (cached.tiles_path/'manifest.json').unlink()
            reload_image = Mock(return_value=image)
            missing = ImageProcessor(None, cached.tiles_path, spatial_transform=transform, image_loader=reload_image)
            self.assertTrue(missing.create_tiles(overwrite=False))
            reload_image.assert_called_once_with()

    def test_image_dimensions_are_read_without_decoding(self):
        import tifffile
        from hd_viewer.startup.data.VisiumHD import VisiumHD
        with tempfile.TemporaryDirectory() as root:
            for suffix, channels in (('png',3), ('png',4), ('tif',3), ('btf',3)):
                with self.subTest(suffix=suffix, channels=channels):
                    path = Path(root)/f'image{channels}.{suffix}'
                    image = np.zeros((7,11,channels), np.uint8)
                    if suffix=='png':
                        Image.fromarray(image).save(path)
                    else:
                        tifffile.imwrite(path, image, bigtiff=suffix=='btf')
                    data = VisiumHD(img_path=str(path), adata_path='unused', positions_path='unused', bins_size=2)
                    with patch('cv2.imread', side_effect=AssertionError('Header lookup must not decode')), \
                         patch('tifffile.TiffPage.asarray', side_effect=AssertionError('Header lookup must not decode')):
                        self.assertEqual(data.image_shape, (7,11,3))
                        self.assertIsNone(data.img)

    def test_tile_source_exposes_the_highest_generated_zoom(self):
        with tempfile.TemporaryDirectory() as root:
            for size in (8, 130):
                processor = ImageProcessor(np.zeros((3, size, 3), np.uint8), root, tile_size=16)
                processor.deploy_server = lambda port, host=None: register_tiles(root)
                source = ImageComponent(processor).generate_tile_source()
                self.assertEqual(source.max_zoom, processor.max_zoom)


if __name__ == "__main__":
    unittest.main()
