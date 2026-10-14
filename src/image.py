
import numpy as np
from PIL import Image

from flask import Flask, send_file, abort
from shapely.affinity import scale as shapely_scale
import os
import math

import threading

class ImageProcessor:
    
    """
    Class to handle image tiles.
    """

    tiles_path = './tiles'
    filename_pattern = 'z{zoom}_x{x}_y{y}.png'

    def __init__(self, image, tiles_path, tile_size=256):
        self.image = np.flipud(image)
        self.tiles_path = tiles_path
        self.tile_size = tile_size

    @property
    def mem(self):
        return self.image.size * self.image.itemsize / 1e9 # GB

    @property
    def max_dim(self):
        return max(self.image.shape[:2])

    @property
    def max_zoom(self):
        """
        Set the maximum zoom level based on the image dimensions.
        """
        return max(0, math.floor(math.log(self.max_dim / self.tile_size, 2)))
    
    @property
    def scale_xfactor(self):
        return self.max_dim / self.image.shape[1]
    
    @property
    def scale_yfactor(self):
        return self.max_dim / self.image.shape[0]
    
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

    def create_tiles(self, overwrite=True, progress_callback=None):
        image = Image.fromarray(self.image)
        base_width, base_height = image.size

        if not os.path.exists(self.tiles_path):
            os.makedirs(self.tiles_path)
        else:
            if len(os.listdir(self.tiles_path)) > 0 and not overwrite:
                return
            self.clean_tiles(rm_dir=False)

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
                if zoom == self.max_zoom:
                    resized_image = image
                else:
                    resized_image = image.resize((scaled_width, scaled_height), Image.LANCZOS)

                for x in range(x_tiles):
                    for y in range(y_tiles):
                        left = x * self.tile_size
                        top = y * self.tile_size
                        right = min(left + self.tile_size, base_width)
                        bottom = min(top + self.tile_size, base_height)

                        tile = resized_image.crop((left, top, right, bottom))

                        if zoom > 0:
                            tile = tile.resize((self.tile_size, self.tile_size), Image.LANCZOS)

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
        if not os.path.exists(self.tiles_path) or overwrite:
            self.create_tiles()

        app = Flask("TileSource")

        @app.route("/<int:zoom>/<int:x>/<int:y>")
        def get_tile(zoom, x, y):
            # Assuming tiles are named based on a fixed format
            filename = self.filename_pattern.format(zoom=zoom, x=x, y=y)
            filepath = os.path.join(self.tiles_path, filename)

            if not os.path.exists(filepath):
                abort(404)

            return send_file(filepath, mimetype='image/png'), {'Access-Control-Allow-Origin':'*'}
        
        def run_app():
            app.run(debug=False, host=host, port=port)

        # Run the Flask app in a separate thread
        threading.Thread(target=run_app).start()
