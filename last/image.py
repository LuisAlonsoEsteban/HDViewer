
from tifffile import imread
from csbdeep.utils import normalize
import numpy as np
from skimage.measure import block_reduce
import hashlib

from . import RAM_LIMIT

class Image:
    def __init__(self, img=None, path:str=None):
        self._img = img
        self.path = path
        self.res_dict = {}

    @property
    def img(self):
        if self._img is None:
            self.read()
        return self._img

    @property
    def width(self):
        return self.img.shape[1]
    
    @property
    def height(self):
        return self.img.shape[0]
    
    @property
    def pixels(self):
        return self.img.shape[0] * self.img.shape[1]
    
    def __getitem__(self, key):
        return self.res_dict[key]
    
    def __setitem__(self, key, value):
        self.res_dict[key] = value

    def clean_cache(self):
        self.res_dict = {}

    def read(self):
        self._img = imread(self.path)

    def normalize(self, min_perc=5, max_perc=95):
        self._img = normalize(self.img, min_perc, max_perc)

    def generate_hash_key(self, bbox):
        return hashlib.md5(str(bbox).encode()).hexdigest()
    
    def get_image_within_memory_limit(self, bbox):
        hash_key = self.generate_hash_key(bbox)
        if hash_key in self.res_dict:
            return self[hash_key]

        x_min, x_max, y_min, y_max = bbox
        cropped_img = self.img[y_min:y_max, x_min:x_max]

        # Calculate the number of bytes per pixel (assuming float32 for now)
        bytes_per_pixel = cropped_img.itemsize
        target_pixel_count = (RAM_LIMIT * (1024**3)) // bytes_per_pixel
    
        # Calculate current pixel count
        current_pixel_count = cropped_img.size

        if current_pixel_count > target_pixel_count:
            # Calculate downsampling factor
            factor = np.sqrt(current_pixel_count / target_pixel_count)
            factor = int(np.ceil(factor))  # Ensure the factor is an integer
            
            print(f"Downsampling by a factor of {factor} to fit within memory limit ({RAM_LIMIT} GB).")

            # Downsample the image
            downsampled_img = block_reduce(cropped_img, block_size=(factor, factor, 1), func=np.mean)
            self[hash_key] = (downsampled_img, factor)
            
            return downsampled_img, factor
        else:
            return cropped_img, None
        
    
class RGBAImage(Image):

    def __init__(self, img=None, path:str=None):
        super().__init__(img, path)
        self._img = None
        self.path = path
        self.res_dict = {}

    def is2d(self):
        return len(self.img.shape) == 2
    
    def is3d(self):
        return len(self.img.shape) == 3
    
    def get_rgba(self, alpha=255):
        if self.is3d() and self.img.shape[2] == 3 and self.img.dtype == np.uint8:
            rgba = np.zeros((self.img.shape[0], self.img.shape[1], 4), dtype=np.uint8)
            rgba[..., :3] = self.img
            rgba[..., 3] = alpha
            return rgba
        else:
            return self.img
    
    def collapse(self):
        rgba = self.get_rgba()
        if rgba.dtype == np.uint8:
            return (rgba[..., 0].astype(np.uint32) << 24) | \
                    (rgba[..., 1].astype(np.uint32) << 16) | \
                    (rgba[..., 2].astype(np.uint32) << 8)  | \
                    (rgba[..., 3].astype(np.uint32))
        else:
            raise ValueError("Image is not in RGBA format or not uint8 type.")
        
    def plot(self, p, bbox, **kwargs):
        rgba = self.collapse()
        cropped_rgba = rgba[bbox[0]:bbox[1], bbox[2]:bbox[3]]
        p.image_rgba(image=[cropped_rgba], x=bbox[0], y=bbox[2], dw=(bbox[3]-bbox[2]), dh=(bbox[1]-bbox[0]), **kwargs)

