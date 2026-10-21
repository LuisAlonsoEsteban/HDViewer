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


if __name__ == "__main__":
    unittest.main()
