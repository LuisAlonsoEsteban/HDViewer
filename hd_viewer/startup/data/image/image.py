import numpy as np
from PIL import Image

from flask import Flask, send_file, abort
import os
import math

import hashlib
import threading
import json
import shutil
import tempfile
import uuid
from pathlib import Path

from .... import settings

_CACHE_LOCK = threading.Lock()

_TILE_SETS = {}                 # key -> (tiles directory, file extension)
_TILE_SERVERS = set()           # ports with a running tile server
_TILE_LOCK = threading.Lock()

TILE_MIMETYPES = {'webp': 'image/webp', 'png': 'image/png'}
TILE_NAME = 'z{zoom}_x{x}_y{y}.{ext}'
WEBP_MAX_DIM = 16383            # libwebp's limit; larger single images fall back to PNG
TILE_CACHE_CONTROL = 'public, max-age=31536000, immutable'


def _read_manifest(path):
    try:
        return (Path(path) / 'manifest.json').read_bytes()
    except OSError:
        return None


def register_tiles(tiles_path):
    """Make a tiles directory servable under a key that changes with its content.

    The key hashes the directory *and* its manifest (format, quality, tile size,
    shape, transform and a per-generation nonce), so every re-tiling gets a new URL
    and tiles can be cached by browsers forever. Older keys of the same directory are
    dropped: their URLs now 404 instead of serving different pixels.
    """
    path = os.path.abspath(str(tiles_path))
    manifest = _read_manifest(path)
    if manifest is None:
        # Unmanaged directory (no manifest): version it by its files' sizes and mtimes.
        entries = sorted((e.name, e.stat().st_size, e.stat().st_mtime_ns) for e in os.scandir(path))
        version, ext = repr(entries).encode(), 'png'
        for name, *_ in entries:
            if name.startswith('z0_x0_y0.'):
                ext = name.rsplit('.', 1)[1]
    else:
        version = manifest
        try:
            ext = json.loads(manifest).get('format', 'png')
        except ValueError:
            ext = 'png'
    key = hashlib.sha1(path.encode() + b'\0' + version).hexdigest()[:16]
    with _TILE_LOCK:
        for old in [k for k, (p, _) in _TILE_SETS.items() if p == path and k != key]:
            del _TILE_SETS[old]
        _TILE_SETS[key] = (path, ext)
    return key


def forget_tiles(tiles_path):
    """Stop serving every key of a directory whose tiles are being replaced."""
    path = os.path.abspath(str(tiles_path))
    with _TILE_LOCK:
        for old in [k for k, (p, _) in _TILE_SETS.items() if p == path]:
            del _TILE_SETS[old]


def tile_extension(key):
    """File extension (and URL suffix) of the tiles registered under `key`."""
    with _TILE_LOCK:
        return _TILE_SETS[key][1]


def make_tile_app():
    app = Flask("TileSource")
    cors = {'Access-Control-Allow-Origin': '*'}

    @app.errorhandler(404)
    def not_found(error):
        # Never cache misses: the same URL may become valid (or the key re-registered).
        return 'Tile not found', 404, {**cors, 'Cache-Control': 'no-store'}

    @app.route("/<key>/<int:zoom>/<int:x>/<int:y>")
    @app.route("/<key>/<int:zoom>/<int:x>/<int:y>.<ext>")
    def get_tile(key, zoom, x, y, ext=None):
        from ....core.tiles import resolve_tile
        found = resolve_tile(key, zoom, x, y, ext)
        if found is None:
            abort(404)
        filepath, mimetype, etag = found
        # The key is content-versioned, so the response for this URL never changes.
        response = send_file(filepath, mimetype=mimetype, etag=etag, conditional=True, max_age=31536000)
        response.headers['Cache-Control'] = TILE_CACHE_CONTROL
        response.headers.update(cors)
        return response

    return app


def ensure_tile_server(host='0.0.0.0', port=4000):
    """Start the process' tile server on `port` once."""
    with _TILE_LOCK:
        if port in _TILE_SERVERS:
            return
        _TILE_SERVERS.add(port)
    app = make_tile_app()
    threading.Thread(target=lambda: app.run(debug=False, host=host, port=port), daemon=True).start()


class ImageProcessor:
    
    """
    Class to handle image tiles.
    """

    tiles_path = './tiles'

    def __init__(self, image, tiles_path, tile_size=256, spatial_transform=None, image_loader=None,
                 tile_format=None, tile_quality=None):
        self.spatial_transform = spatial_transform
        self._source_image = image
        self._image_loader = image_loader
        self._image = None
        self._single_image = False
        self.tiles_path = Path(tiles_path)
        self.tile_size = tile_size
        self.tile_format = (tile_format or settings.TILE_FORMAT).lower()
        if self.tile_format not in TILE_MIMETYPES:
            raise ValueError(f"Unsupported tile format {self.tile_format!r}; use one of {sorted(TILE_MIMETYPES)}")
        self.tile_quality = int(settings.TILE_QUALITY if tile_quality is None else tile_quality)

    @property
    def tile_extension(self):
        """Format actually written: WebP cannot hold a single image above 16383 px."""
        if self.tile_format == 'webp' and self._single_image and self.max_dim > WEBP_MAX_DIM:
            return 'png'
        return self.tile_format

    @property
    def filename_pattern(self):
        return TILE_NAME.replace('{ext}', self.tile_extension)

    def _save_tile(self, tile, path):
        """Encode one tile; alpha is kept only where the padding makes it transparent."""
        if tile.mode not in ('RGB', 'RGBA'):
            tile = tile.convert('RGBA' if tile.has_transparency_data else 'RGB')
        if tile.mode == 'RGBA' and tile.getchannel('A').getextrema() == (255, 255):
            tile = tile.convert('RGB')
        if self.tile_extension == 'webp':
            tile.save(path, 'WEBP', quality=self.tile_quality, method=4)
        else:
            tile.save(path, 'PNG')

    @property
    def image(self):
        """Decode and pad pixels only when a cache miss actually needs them."""
        if self._image is None:
            source = self._source_image
            if source is None:
                if self._image_loader is None:
                    raise ValueError('Image pixels or an image loader are required.')
                source = self._image_loader()
            if source is None:
                raise ValueError('Cannot read the tissue image.')
            self._image = (self.spatial_transform.warp_image(source) if self.spatial_transform
                           else np.flipud(source))
            self._source_image = None
        return self._image

    @property
    def shape(self):
        if self.spatial_transform is not None:
            size = self.spatial_transform.output_size
            return (size, size, 4)  # pad_image/warp_image produces uint8 RGBA
        return self.image.shape

    @property
    def mem(self):
        if self.spatial_transform is not None:
            return int(np.prod(self.shape)) / 1e9
        return self.image.size * self.image.itemsize / 1e9

    @property
    def max_dim(self):
        return max(self.shape[:2])

    @property
    def max_zoom(self):
        """
        Set the maximum zoom level based on the image dimensions.
        """
        return max(0, math.floor(math.log(self.max_dim / self.tile_size, 2)))
    
    @property
    def scale_xfactor(self):
        return self.max_dim / self.shape[1]
    
    @property
    def scale_yfactor(self):
        return self.max_dim / self.shape[0]
    
    def clean_tiles(self, rm_dir=True):
        """
        Clean the tiles directory by removing all files.
        """
        if os.path.exists(self.tiles_path):
            for file_name in os.listdir(self.tiles_path):
                file_path = os.path.join(self.tiles_path, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            if rm_dir:
                os.rmdir(self.tiles_path)

    def create_tiles(self, overwrite=True, progress_callback=None, single_image=False):
        self._single_image = single_image
        ext = self.tile_extension
        # Format and quality are part of the cache identity: PNG caches built by older
        # versions (version 2, no format) are regenerated rather than served as WebP.
        spec = dict(version=3, shape=list(self.shape), tile_size=self.tile_size,
                    transform=getattr(self.spatial_transform, "fingerprint", None), single_image=single_image,
                    format=ext, quality=self.tile_quality if ext == 'webp' else None)
        self.tiles_path.parent.mkdir(parents=True, exist_ok=True)

        with _CACHE_LOCK:
            marker = self.tiles_path / "manifest.json"
            if not overwrite and marker.exists():
                try:
                    cached = json.loads(marker.read_text())
                    # Reuse any cache with the same pixels, whatever its format or
                    # quality (e.g. PNG caches from older versions): re-tiling a large
                    # slide costs far more than serving its existing tiles.
                    geometry = ("shape", "tile_size", "transform", "single_image")
                    if all(cached.get(name) == spec[name] for name in geometry):
                        return False
                except (ValueError, OSError):
                    pass

            temporary = Path(tempfile.mkdtemp(prefix=".tiles-", dir=self.tiles_path.parent))
            original = self.tiles_path

            try:
                self.tiles_path = temporary
                if single_image:
                    self._save_tile(Image.fromarray(self.image),
                                    temporary / self.filename_pattern.format(zoom=0, x=0, y=0))
                    if progress_callback:
                        progress_callback(1, 1)
                else:
                    self._write_tiles(progress_callback)

                # A fresh nonce per generation gives re-tiled content a new URL key.
                (temporary / "manifest.json").write_text(
                    json.dumps(dict(spec, generation=uuid.uuid4().hex), sort_keys=True))
                forget_tiles(original)
                
                backup = None
                if original.exists():
                    backup = Path(tempfile.mkdtemp(prefix=".old-tiles-", dir=original.parent))
                    backup.rmdir()
                    original.rename(backup)
                
                try:
                    temporary.rename(original)
                except Exception:
                    if backup is not None:
                        backup.rename(original)
                    raise

                if backup is not None:
                    shutil.rmtree(backup)
                    
            finally:
                self.tiles_path = original
                if temporary.exists():
                    shutil.rmtree(temporary)
        return True

    def _write_tiles(self, progress_callback=None):
        image = Image.fromarray(self.image)

        total_steps = 0
        zoom_tile_counts = []

        for zoom in range(self.max_zoom + 1):
            factor = 2 ** zoom
            scaled_width = self.tile_size * factor
            scaled_height = self.tile_size * factor
            x_tiles = math.ceil(scaled_width / self.tile_size)
            y_tiles = math.ceil(scaled_height / self.tile_size)
            count = x_tiles * y_tiles
            zoom_tile_counts.append((zoom, x_tiles, y_tiles, scaled_width, scaled_height))
            total_steps += count

        done = 0

        try:
            for zoom, x_tiles, y_tiles, scaled_width, scaled_height in zoom_tile_counts:
                if image.size == (scaled_width, scaled_height):
                    resized_image = image
                else:
                    resized_image = image.resize((scaled_width, scaled_height), Image.LANCZOS)

                for x in range(x_tiles):
                    for y in range(y_tiles):
                        left = x * self.tile_size
                        top = y * self.tile_size
                        right = left + self.tile_size
                        bottom = top + self.tile_size

                        tile = resized_image.crop((left, top, right, bottom))

                        tile_path = os.path.join(
                            self.tiles_path,
                            self.filename_pattern.format(zoom=zoom, x=x, y=y)
                        )
                        self._save_tile(tile, tile_path)

                        done += 1
                        if progress_callback is not None:
                            progress_callback(done, total_steps)
        
        except Exception:
            self.clean_tiles()
            raise


    def deploy_server(self, host='0.0.0.0', port=4000, overwrite=False):
        """Serve these tiles; returns the key of their URL path (/<key>/<z>/<x>/<y>.<ext>).

        One tile server per process and port serves every image: starting one per
        dataset failed silently once the port was taken, and the browser was then
        shown the first dataset's image for every later one.
        """
        
        self.create_tiles(overwrite=overwrite, single_image=self._single_image)
        self.tile_key = register_tiles(self.tiles_path)
        ensure_tile_server(host, port)
        return self.tile_key
