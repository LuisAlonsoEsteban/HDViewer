
from PIL import Image

from flask import Flask, send_file, abort

import os
from tqdm import tqdm
import math

import threading

PORT = 5001
FLASK_ROUTE = "/<int:zoom>/<int:x>/<int:y>"

def get_url():
    f"http://localhost:{PORT}"

def get_arguments():
    return "/".join([d.split(":")[-1].capitalize() for d in FLASK_ROUTE.split("/")])

def get_request():
    url = get_url()
    arguments = get_arguments()
    return f"{url}{arguments}"


class ImageTiles:
    
    """
    Class to handle image tiles.
    """

    tiles_path = './'
    filename_pattern = 'z{zoom}_x{x}_y{y}.png'
    
    PORT = 4000
    QUERY = "/<int:zoom>/<int:x>/<int:y>"
    URL = f"http://localhost:{PORT}{QUERY}"

    def __init__(self, image, tiles_path, tile_size=256):
        self.image = image
        self.tiles_path = tiles_path
        self.tile_size = tile_size

    @property
    def max_zoom(self):
        """
        Set the maximum zoom level based on the image dimensions.
        """
        base_width, base_height = self.image.shape[:2]
        max_dim = max(base_width, base_height)
        return max(0, math.floor(math.log(max_dim / self.tile_size, 2)))

    def clean_tiles(self):
        """
        Clean the tiles directory by removing all files.
        """
        if os.path.exists(self.tiles_path):
            for file_name in os.listdir(self.tiles_path):
                file_path = os.path.join(self.tiles_path, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            os.rmdir(self.tiles_path)

    def create_tiles(self):
        image = Image.fromarray(self.image)
        base_width, base_height = image.size

        if not os.path.exists(self.tiles_path):
            os.makedirs(self.tiles_path)
        else:
            self.clean_tiles()

        try:
            for zoom in range(self.max_zoom + 1):
                factor = 2 ** zoom
                scaled_width = self.tile_size * factor
                scaled_height = self.tile_size * factor

                if zoom == self.max_zoom:
                    resized_image = image
                else:
                    resized_image = image.resize((scaled_width, scaled_height), Image.LANCZOS)
                
                # Determine the number of tiles at this zoom level
                x_tiles = math.ceil(scaled_width / self.tile_size)
                y_tiles = math.ceil(scaled_height / self.tile_size)

                for x in tqdm(range(x_tiles), desc=f"Processing zoom level {zoom}, out of {self.max_zoom} zoom levels"):
                    for y in range(y_tiles):
                        left = x * self.tile_size
                        top = y * self.tile_size
                        right = min(left + self.tile_size, base_width)
                        bottom = min(top + self.tile_size, base_height)

                        # Crop directly from the base image
                        tile = resized_image.crop((left, top, right, bottom))
                        
                        # Resize only if necessary (for zoom levels > 0)
                        if zoom > 0:
                            tile = tile.resize((self.tile_size, self.tile_size), Image.LANCZOS)
                        
                        # Save the tile
                        tile_path = os.path.join(self.tiles_path, self.filename_pattern.format(zoom=zoom, x=x, y=y))
                        tile.save(tile_path)
        except:
            self.clean_tiles()
            raise KeyError


    def deploy_server(self, host='0.0.0.0', overwrite=False):
        if not os.path.exists(self.tiles_path) or overwrite:
            self.create_tiles()

        app = Flask("TileSource")

        @app.route(FLASK_ROUTE)
        def get_tile(zoom, x, y):
            # Assuming tiles are named based on a fixed format
            filename = self.filename_pattern.format(zoom=zoom, x=x, y=y)
            filepath = os.path.join(self.tiles_path, filename)

            if not os.path.exists(filepath):
                abort(404)

            return send_file(filepath, mimetype='image/png'), {'Access-Control-Allow-Origin':'*'}
        
        def run_app():
            app.run(debug=False, host=host, port=PORT)

        # Run the Flask app in a separate thread
        threading.Thread(target=run_app).start()

        return get_request()