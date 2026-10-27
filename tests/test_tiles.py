import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image

from bokeh_server.components.image import ImageComponent
from src.image import ImageProcessor, make_tile_app, register_tiles


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
            processor.deploy_server = lambda port: register_tiles(processor.tiles_path)   # no real server
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

    def test_tile_source_exposes_the_highest_generated_zoom(self):
        with tempfile.TemporaryDirectory() as root:
            for size in (8, 130):
                processor = ImageProcessor(np.zeros((3, size, 3), np.uint8), root, tile_size=16)
                processor.deploy_server = lambda port: register_tiles(root)
                source = ImageComponent(processor).generate_tile_source()
                self.assertEqual(source.max_zoom, processor.max_zoom)


if __name__ == "__main__":
    unittest.main()
