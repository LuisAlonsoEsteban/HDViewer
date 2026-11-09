import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
from PIL import Image

from hd_viewer.api import figures
from hd_viewer.core.tiles import resolve_tile
from hd_viewer.startup.data.image.image import ImageProcessor, register_tiles, tile_extension


def tile_client():
    """The API's tile route (GET /api/tiles/{key}/{zoom}/{x}/{y}[.ext]) on its own."""
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from hd_viewer.api.routers import tiles
    app = FastAPI()
    app.include_router(tiles.router, prefix="/api")
    return TestClient(app)


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
            client = tile_client()
            pixel = lambda key: np.asarray(Image.open(
                io.BytesIO(client.get(f"/api/tiles/{key}/0/0/0").content)))[0, 0, 0]
            self.assertEqual(pixel(brain), 10)
            self.assertEqual(pixel(ovary), 200)      # not the first dataset's image
            self.assertEqual(client.get("/api/tiles/unknown/0/0/0").status_code, 404)
            self.assertEqual(client.get(f"/api/tiles/{ovary}/3/0/0").status_code, 404)
            self.assertEqual(client.get(f"/api/tiles/{ovary}/0/0/0.png").status_code, 200)
            self.assertEqual(client.get(f"/api/tiles/{ovary}/0/0/0.webp").status_code, 404)   # wrong extension
            self.assertEqual(client.get(f"/api/tiles/{ovary}/0/0/x.png").status_code, 404)    # malformed
            path, mimetype, etag = resolve_tile(ovary, 0, 0, 0)
            self.assertEqual((Path(path).parent, mimetype), (Path(root) / "ovary", "image/png"))
            self.assertIn(ovary, etag)
            self.assertIsNone(resolve_tile(ovary, 0, 0, 0, "webp"))

    def test_image_urls_carry_the_key(self):
        with tempfile.TemporaryDirectory() as root:
            processor = ImageProcessor(np.zeros((8, 8, 3), np.uint8), tile_dir(root, "img", 0))
            key = register_tiles(processor.tiles_path)
            self.assertEqual(figures.tile_template(processor, "/app"),
                             f"/app/api/tiles/{key}/{{Z}}/{{X}}/{{Y}}.png")
            plot = SimpleNamespace(calls=[], image_url=lambda **kw: plot.calls.append(kw))
            figures._image(plot, processor, True, "")
            self.assertEqual(plot.calls[0]["url"], [f"/api/tiles/{key}/0/0/0.png"])
            self.assertEqual(tile_client().get(plot.calls[0]["url"][0]).status_code, 200)


class TileGenerationTests(unittest.TestCase):
    def test_every_zoom_covers_rectangular_and_small_images(self):
        for height, width in ((3, 130), (130, 3), (71, 101), (5, 7)):
            with self.subTest(shape=(height, width)), tempfile.TemporaryDirectory() as root:
                image = np.random.default_rng(0).integers(0, 256, (height, width, 3), dtype=np.uint8)
                processor = ImageProcessor(image, Path(root) / "tiles", tile_size=16, tile_format="png")
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
                processor = ImageProcessor(np.zeros((3, size, 3), np.uint8), Path(root) / str(size), tile_size=16)
                processor.create_tiles()
                from bokeh.plotting import figure
                plot = figure()
                figures._image(plot, processor, False, "")
                source = plot.select_one({"name": "image:tiles"}).tile_source
                self.assertEqual(source.max_zoom, processor.max_zoom)
                self.assertTrue(source.url.endswith(f"/{register_tiles(processor.tiles_path)}/{{Z}}/{{X}}/{{Y}}.webp"),
                                source.url)


def padded(height, width, size):
    """An opaque gradient in the corner of a transparent square, like SpatialTransform.warp_image."""
    canvas = np.zeros((size, size, 4), np.uint8)
    yy, xx = np.mgrid[:height, :width]
    canvas[:height, :width, :3] = np.stack([4 * xx, 4 * yy, 2 * (xx + yy)], -1).clip(0, 255)
    canvas[:height, :width, 3] = 255
    return canvas


class TileFormatTests(unittest.TestCase):
    def test_webp_tiles_keep_alpha_only_where_padding_is_transparent(self):
        with tempfile.TemporaryDirectory() as root:
            # np.flipud puts the opaque rows at the bottom: tiles y=1 are opaque, y=0 transparent.
            processor = ImageProcessor(padded(32, 64, 64), Path(root) / "t", tile_size=32,
                                       tile_format="webp", tile_quality=80)
            processor.create_tiles()
            self.assertEqual(processor.max_zoom, 1)
            def tile(x, y):
                with Image.open(processor.tiles_path / f"z1_x{x}_y{y}.webp") as image:
                    image.load()
                    return image
            self.assertEqual(tile(0, 1).format, "WEBP")
            self.assertEqual(tile(0, 1).mode, "RGB")           # fully opaque -> no alpha channel
            self.assertEqual(tile(0, 0).mode, "RGBA")          # padding -> alpha kept
            self.assertEqual(tile(0, 0).getchannel("A").getextrema(), (0, 0))
            expected = np.flipud(padded(32, 64, 64))[32:, :32, :3].astype(int)
            self.assertLess(np.abs(np.asarray(tile(0, 1)).astype(int) - expected).mean(), 4)   # lossy but close
            manifest = json.loads((processor.tiles_path / "manifest.json").read_text())
            self.assertEqual((manifest["format"], manifest["quality"]), ("webp", 80))

    def test_server_sends_format_mimetype_and_immutable_cache_headers(self):
        with tempfile.TemporaryDirectory() as root:
            client = tile_client()
            for fmt in ("webp", "png"):
                with self.subTest(fmt=fmt):
                    processor = ImageProcessor(padded(20, 20, 40), Path(root) / fmt, tile_size=16, tile_format=fmt)
                    processor.create_tiles()
                    key = register_tiles(processor.tiles_path)
                    self.assertEqual(tile_extension(key), fmt)
                    response = client.get(f"/api/tiles/{key}/1/0/0.{fmt}")
                    self.assertEqual(response.status_code, 200)
                    self.assertEqual(response.headers["Content-Type"], f"image/{fmt}")
                    with Image.open(io.BytesIO(response.content)) as image:
                        self.assertEqual(image.format, fmt.upper())
                    self.assertEqual(response.headers["Cache-Control"], "public, max-age=31536000, immutable")
                    self.assertEqual(response.headers["Access-Control-Allow-Origin"], "*")
                    etag = response.headers["ETag"]
                    self.assertIn(key, etag)
                    again = client.get(f"/api/tiles/{key}/1/0/0.{fmt}", headers={"If-None-Match": etag})
                    self.assertEqual(again.status_code, 304)
                    self.assertIn("immutable", again.headers["Cache-Control"])
                    for missing in (f"/api/tiles/{key}/9/0/0.{fmt}", f"/api/tiles/nokey/0/0/0.{fmt}"):
                        miss = client.get(missing)
                        self.assertEqual(miss.status_code, 404)
                        self.assertNotIn("immutable", miss.headers.get("Cache-Control", ""))
                        self.assertIn("no-store", miss.headers["Cache-Control"])

    def test_key_changes_with_content_format_and_regeneration(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "tiles"
            image = padded(20, 20, 40)
            ImageProcessor(image, path, tile_size=16, tile_format="webp").create_tiles()
            first = register_tiles(path)
            self.assertEqual(register_tiles(path), first)               # stable while unchanged
            # Cache hit (same spec, overwrite off) keeps the key and the files.
            self.assertFalse(ImageProcessor(image, path, tile_size=16, tile_format="webp").create_tiles(overwrite=False))
            self.assertEqual(register_tiles(path), first)
            keys = {first}
            # Format/quality alone: the existing tiles are reused; a rebuild needs overwrite
            for kwargs in (dict(tile_format="webp", tile_quality=50), dict(tile_format="png")):
                with self.subTest(**kwargs):
                    processor = ImageProcessor(image, path, **{"tile_size": 16, **kwargs})
                    self.assertFalse(processor.create_tiles(overwrite=False))
                    self.assertTrue(processor.create_tiles(overwrite=True))
                    key = register_tiles(path)
                    self.assertNotIn(key, keys)
                    keys.add(key)
            # Different pixels (tile size): rebuilt even without overwrite
            processor = ImageProcessor(image, path, tile_size=8, tile_format="png")
            self.assertTrue(processor.create_tiles(overwrite=False))
            self.assertNotIn(register_tiles(path), keys)
            # Same settings, overwrite=True (e.g. a new image at the same path): new key, old one stops serving.
            before = register_tiles(path)
            ImageProcessor(image[::-1].copy(), path, tile_size=8, tile_format="png").create_tiles(overwrite=True)
            client = tile_client()
            self.assertEqual(client.get(f"/api/tiles/{before}/0/0/0.png").status_code, 404)
            after = register_tiles(path)
            self.assertNotEqual(after, before)
            self.assertEqual(client.get(f"/api/tiles/{after}/0/0/0.png").status_code, 200)

    def test_legacy_png_cache_is_reused_without_re_tiling(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "tiles"
            image = padded(20, 20, 40)
            legacy = ImageProcessor(image, path, tile_size=16, tile_format="png")
            legacy.create_tiles()
            manifest = json.loads((path / "manifest.json").read_text())
            for field in ("format", "quality", "generation"):
                manifest.pop(field)
            (path / "manifest.json").write_text(json.dumps(dict(manifest, version=2)))   # pre-WebP manifest
            self.assertEqual(tile_extension(register_tiles(path)), "png")    # still served correctly as PNG
            webp = ImageProcessor(image, path, tile_size=16, tile_format="webp")
            with patch.object(ImageProcessor, "_write_tiles") as write:
                self.assertFalse(webp.create_tiles(overwrite=False))     # same pixels: reused
                write.assert_not_called()
            self.assertEqual(tile_extension(register_tiles(path)), "png")
            self.assertTrue(webp.create_tiles(overwrite=True))           # explicit rebuild → WebP
            self.assertTrue((path / "z0_x0_y0.webp").exists())
            self.assertFalse((path / "z0_x0_y0.png").exists())
            self.assertEqual(tile_extension(register_tiles(path)), "webp")

    def test_single_image_uses_webp_and_falls_back_to_png_above_webp_limit(self):
        with tempfile.TemporaryDirectory() as root:
            processor = ImageProcessor(padded(20, 30, 40), Path(root) / "small", tile_size=16, tile_format="webp")
            processor.create_tiles(single_image=True)
            with Image.open(processor.tiles_path / "z0_x0_y0.webp") as image:
                self.assertEqual((image.format, image.size, image.mode), ("WEBP", (40, 40), "RGBA"))
            huge = ImageProcessor(None, Path(root) / "huge", tile_format="webp")
            huge._single_image = True
            with patch.object(ImageProcessor, "shape", (20000, 20000, 4)):
                self.assertEqual(huge.tile_extension, "png")

    def test_unknown_format_is_rejected(self):
        with self.assertRaises(ValueError):
            ImageProcessor(None, "unused", tile_format="gif")


if __name__ == "__main__":
    unittest.main()
