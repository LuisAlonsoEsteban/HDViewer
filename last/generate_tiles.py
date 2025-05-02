
from tifffile import imread
from PIL import Image
import numpy as np

import os
from tqdm import tqdm
import math

def make_square(img:np.ndarray):
    height, width = img.shape[:2]
    if width > height:
        padding = (width - height) // 2
        padded_array = np.pad(img, ((padding, width - height - padding), (0, 0), (0, 0)), mode='edge')
    elif height > width:
        padding = (height - width) // 2
        padded_array = np.pad(img, ((0, 0), (padding, height - width - padding), (0, 0)), mode='edge')
    else:
        padded_array = img

    return padded_array


def create_tiles(image_path, tile_directory="tiles", tile_size=(256, 256)):
    image = imread(image_path)
    # image = make_square(image)

    image = Image.fromarray(image)
    base_width, base_height = image.size

    max_dim = max(base_width, base_height)
    max_zoom = math.floor(math.log(max_dim / tile_size[0], 2))
    max_zoom = max(0, max_zoom)

    tile_directory = os.path.join(os.path.dirname(image_path), tile_directory)
    
    if not os.path.exists(tile_directory):
        os.makedirs(tile_directory)
    else:
        for file_name in os.listdir(tile_directory):
            file_path = os.path.join(tile_directory, file_name)
            if os.path.isfile(file_path):
                os.remove(file_path)

    for zoom in range(max_zoom + 1):
        factor = 2 ** zoom
        scaled_width = tile_size[0] * factor
        scaled_height = tile_size[1] * factor

        if zoom == max_zoom:
            resized_image = image
        else:
            resized_image = image.resize((scaled_width, scaled_height), Image.LANCZOS)
        
        # Determine the number of tiles at this zoom level
        x_tiles = math.ceil(scaled_width / tile_size[0])
        y_tiles = math.ceil(scaled_height / tile_size[1])

        for x in tqdm(range(x_tiles), desc=f"Processing zoom level {zoom}, out of {max_zoom} zoom levels"):
            for y in range(y_tiles):
                left = x * tile_size[0]
                top = y * tile_size[1]
                right = min(left + tile_size[0], base_width)
                bottom = min(top + tile_size[1], base_height)

                # Crop directly from the base image
                tile = resized_image.crop((left, top, right, bottom))
                
                # Resize only if necessary (for zoom levels > 0)
                if zoom > 0:
                    tile = tile.resize(tile_size, Image.LANCZOS)
                
                # Save the tile
                tile_path = f"{tile_directory}/z{zoom}_x{x}_y{y}.png"
                tile.save(tile_path)
    
if __name__ == "__main__":
    # Example usage
    create_tiles('/scratch/lalonsoeste/PhD/SpatialTranscriptomics/data/spatial/10x/VisiumHD/Mouse intestine/Visium_HD_Mouse_Small_Intestine_tissue_image.btf')
