
import numpy as np
from PIL import Image

from flask import Flask, send_file, abort
import os
import math

import threading
import json
import shutil
import tempfile
from pathlib import Path

_CACHE_LOCK = threading.Lock()

_TILE_SETS = {}                 # key -> tiles directory
_TILE_SERVERS = set()           # ports with a running tile server
_TILE_LOCK = threading.Lock()


def register_tiles(tiles_path):
    """Make a tiles directory servable; its key changes with the directory."""
    import hashlib
    path = os.path.abspath(str(tiles_path))
    key = hashlib.sha1(path.encode()).hexdigest()[:12]
    with _TILE_LOCK:
        _TILE_SETS[key] = path
    return key


def make_tile_app():
    app = Flask("TileSource")

    @app.route("/<key>/<int:zoom>/<int:x>/<int:y>")
    def get_tile(key, zoom, x, y):
        tiles_path = _TILE_SETS.get(key)
        if tiles_path is None:
            abort(404)
        filepath = os.path.join(tiles_path, ImageProcessor.filename_pattern.format(zoom=zoom, x=x, y=y))
        if not os.path.exists(filepath):
            abort(404)
        return send_file(filepath, mimetype='image/png'), {'Access-Control-Allow-Origin': '*'}

    return app


def ensure_tile_server(host='0.0.0.0', port=5001):
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
    filename_pattern = 'z{zoom}_x{x}_y{y}.png'

    def __init__(self, image, tiles_path, tile_size=256, spatial_transform=None, image_loader=None):
        self.spatial_transform = spatial_transform
        self._source_image = image
        self._image_loader = image_loader
        self._image = None
        self._single_image = False
        self.tiles_path = Path(tiles_path)
        self.tile_size = tile_size

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
        spec = dict(version=2, shape=list(self.shape), tile_size=self.tile_size,
                    transform=getattr(self.spatial_transform, "fingerprint", None), single_image=single_image)
        self.tiles_path.parent.mkdir(parents=True, exist_ok=True)

        with _CACHE_LOCK:
            marker = self.tiles_path / "manifest.json"
            if not overwrite and marker.exists():
                try:
                    if json.loads(marker.read_text()) == spec:
                        return False
                except (ValueError, OSError):
                    pass

            temporary = Path(tempfile.mkdtemp(prefix=".tiles-", dir=self.tiles_path.parent))
            original = self.tiles_path

            try:
                self.tiles_path = temporary
                if single_image:
                    Image.fromarray(self.image).save(temporary / self.filename_pattern.format(zoom=0, x=0, y=0))
                    if progress_callback:
                        progress_callback(1, 1)
                else:
                    self._write_tiles(progress_callback)

                (temporary / "manifest.json").write_text(json.dumps(spec, sort_keys=True))
                
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
                        tile.save(tile_path)

                        done += 1
                        if progress_callback is not None:
                            progress_callback(done, total_steps)
        
        except Exception:
            self.clean_tiles()
            raise


    def deploy_server(self, host='0.0.0.0', port=5001, overwrite=False):
        """Serve these tiles; returns the key of their URL path (/<key>/<z>/<x>/<y>).

        One tile server per process and port serves every image: starting one per
        dataset failed silently once the port was taken, and the browser was then
        shown the first dataset's image for every later one.
        """
        
        self.create_tiles(overwrite=overwrite, single_image=self._single_image)
        self.tile_key = register_tiles(self.tiles_path)
        ensure_tile_server(host, port)
        return self.tile_key
