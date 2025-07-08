
import numpy as np
from PIL import Image

from flask import Flask, send_file, abort

import os
from tqdm import tqdm
import math

import threading
import pickle

from tqdm import tqdm
import geopandas as gpd
from shapely import Polygon

class ImageTiles:
    
    """
    Class to handle image tiles.
    """

    tiles_path = './tiles'
    filename_pattern = 'z{zoom}_x{x}_y{y}.png'

    segmentation_path = './StarDist'

    def __init__(self, image, tiles_path, segm_path, tile_size=256):
        self.image = np.flipud(image)
        self.tiles_path = tiles_path
        self.tile_size = tile_size
        self.segm_path = segm_path

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

    def segment_image(self, axes='YXC', block_size=4096, prob_thresh=0.01,nms_thresh=0.001, min_overlap=128, context=128, normalizer=None, n_tiles=(4,4,1)):
        labels_path = os.path.join(self.segm_path, 'labels.pickle')
        polys_path = os.path.join(self.segm_path, 'polys.pickle')

        if os.path.exists(labels_path) and os.path.exists(polys_path):
            # Load the labels and polys from files
            with open(labels_path, 'rb') as f:
                labels = pickle.load(f)
            with open(polys_path, 'rb') as f:
                polys = pickle.load(f)
            print("Loaded labels and polys from disk.")
        else:           
            from stardist.models import StarDist2D 
            model = StarDist2D.from_pretrained('2D_versatile_he')
            labels, polys = model.predict_instances_big(self.image, axes=axes, block_size=block_size, prob_thresh=prob_thresh, nms_thresh=nms_thresh, min_overlap=min_overlap, context=context, normalizer=normalizer, n_tiles=n_tiles)
            
            with open(labels_path, 'wb') as f:
                pickle.dump(labels, f, protocol=pickle.HIGHEST_PROTOCOL)
            with open(polys_path, 'wb') as f:
                pickle.dump(polys, f, protocol=pickle.HIGHEST_PROTOCOL)
            print("Computed labels and polys and saved to disk.")

        return labels, polys

    def get_cells_gdf(self, polys):
        geometries = []
        for nuclei in tqdm(range(len(polys['coord']))):
            coords = [(y * self.scale_xfactor, x * self.scale_yfactor) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]
            geometries.append(Polygon(coords))

        gdf = gpd.GeoDataFrame(geometry=geometries)
        gdf['index'] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
        gdf.index = gdf['index']
        gdf.index.name = None
        return gdf
