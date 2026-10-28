
import os
import re

import pandas as pd
import scanpy as sc

import json
from pathlib import Path

class VisiumHD:

    img = None
    _adata = None
    _df_tissue_positions = None

    def __init__(
            self,
            dir_base:str = None, 
            img_file:str = None,
            img_path:str = None,
            adata_path:str = None,
            positions_path:str = None,
            bins_size:str = None
            ):
        
        self.dir_base = dir_base
        self.img_path = img_path or os.path.join(self.dir_base, img_file)

        self.adata_path = adata_path or self.get_default_path(bins_size, data='adata')
        self.positions_path = positions_path or self.get_default_path(bins_size, data='positions')

        if bins_size:
            self.bins_size = bins_size
        else:
            assert self.get_bins_size(adata_path) == self.get_bins_size(positions_path)
            self.bins_size = self.get_bins_size(adata_path)

        self.pixels_per_bin = self.get_pixels_per_bin() 

    def read_scalefactors(self):
        if self.dir_base and self.bins_size:
            filename = os.path.join(self.dir_base, f"binned_outputs/square_{self.bins_size:03}um/spatial",
                                    "scalefactors_json.json")
        elif self.positions_path:
            filename = os.path.join(Path(self.positions_path).resolve().parent, "scalefactors_json.json")
        else:
            return None
        if not os.path.exists(filename):
            return None
        with open(filename, 'r') as f:
            return json.load(f)

    @property
    def image_scale(self):
        """Image pixels per full-resolution pixel.

        Positions are in full-resolution pixels. Space Ranger's down-scaled
        tissue_hires_image / tissue_lowres_image need its scale factor; any other
        image is taken to be the full-resolution one (1).
        """
        name = Path(str(self.img_path)).name.lower()
        scalefactors = self.read_scalefactors() or {}
        if "hires" in name and "tissue_hires_scalef" in scalefactors:
            return float(scalefactors["tissue_hires_scalef"])
        if "lowres" in name and "tissue_lowres_scalef" in scalefactors:
            return float(scalefactors["tissue_lowres_scalef"])
        return 1.0

    def get_pixels_per_bin(self):
        if self.dir_base and self.bins_size:
            filename = os.path.join(
                self.dir_base,
                f"binned_outputs/square_{self.bins_size:03}um/spatial",
                "scalefactors_json.json"
            )
        elif self.positions_path:
            filename = os.path.join(
                Path(self.positions_path).resolve().parent,
                "scalefactors_json.json"
            )
        else:
            return
        
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                scalefactors = json.load(f)
            return scalefactors["spot_diameter_fullres"]
        else:
            return

    @staticmethod
    def get_bins_size(base_dir:str) -> int:
        try:
            bins_size_comp = base_dir.split(os.sep)[-2]
            return int(re.search(r'\d{3}', bins_size_comp).group())
        except IndexError:
            raise ValueError(f"Path {base_dir} has no parent directory on which bins size can be retrieved.")
        except AttributeError:
            raise ValueError(f"Path {base_dir} has no bins size in it. You must provide it.")
        
    @staticmethod
    def get_format(path:str) -> str:
        return Path(path).suffix[1:]
        
    def get_available_sizes(self):
        if self.dir_base and os.path.exists(self.dir_base):
            base_path = os.path.join(self.dir_base, "binned_outputs")
            pattern = re.compile(r'square_(\d{3})um')
            available_sizes = []
            for entry in os.scandir(base_path):
                if entry.is_dir() and pattern.search(entry.name):
                    size = pattern.search(entry.name).group(1)
                    available_sizes.append(int(size))
            return available_sizes
        else:
            raise ValueError

    def get_default_path(self, bins_size:str, data='adata'):
        base_path = os.path.join(self.dir_base, f"binned_outputs/square_{bins_size:03}um")
        if data == 'adata':
            return os.path.join(base_path, "filtered_feature_bc_matrix.h5")
        elif data == 'positions':
            return os.path.join(base_path, "spatial", "tissue_positions.parquet")
        elif data == 'both':
            return self.get_default_path(bins_size, data='adata'), self.get_default_path(bins_size, data='positions')
        
    @property
    def image_shape(self):
        if self.img is not None:
            return self.img.shape

        image_format = VisiumHD.get_format(str(self.img_path)).lower()

        if image_format in ('btf', 'tif', 'tiff'):
            from tifffile import TiffFile
            with TiffFile(self.img_path) as image:
                if image_format == 'btf':
                    return image.series[0].shape
                # OpenCV reads the first TIFF page in colour by default.
                page = image.pages[0]
                return (page.imagelength, page.imagewidth, 3)

        if image_format == 'png':
            import struct
            with open(self.img_path, 'rb') as image:
                header = image.read(24)
            if len(header) != 24 or header[:8] != b'\x89PNG\r\n\x1a\n' or header[12:16] != b'IHDR':
                raise ValueError('Invalid PNG image header.')
            width, height = struct.unpack('>II', header[16:24])
            return (height, width, 3)

        raise ValueError('Image format not supported')

    def read_img(self):
        self.img_format = VisiumHD.get_format(str(self.img_path)).lower()
        if self.img_format == "btf":
            from tifffile import imread
        elif self.img_format in ("png", "tif", "tiff"):
            from cv2 import imread
        else:
            raise ValueError("Image format not supported")
        
        self.img = imread(self.img_path)

    def read_positions(self):
        self._df_tissue_positions = pd.read_parquet(self.positions_path)
        if 'barcode' in self._df_tissue_positions.columns:
            self._df_tissue_positions = self._df_tissue_positions.set_index('barcode')
        else:
            raise ValueError("No valid index column found in tissue positions DataFrame")
        self._df_tissue_positions['index'] = self._df_tissue_positions.index
        self.df_tissue_positions = self._df_tissue_positions

    def read_adata(self):
        self.adata = sc.read_10x_h5(self.adata_path)
        self.adata.var_names_make_unique()

    @property
    def adata(self):
        return self._adata

    @adata.setter
    def adata(self, value):
        self._adata = value
        if self.is_data_loaded():
            self.on_data_loaded()

    @property
    def df_tissue_positions(self):
        return self._df_tissue_positions

    @df_tissue_positions.setter
    def df_tissue_positions(self, value):
        self._df_tissue_positions = value
        if self.is_data_loaded():
            self.on_data_loaded()

    def is_data_loaded(self):
        return self._adata is not None and self._df_tissue_positions is not None

    def merge_data(self):
        self._adata.obs =  pd.merge(self._adata.obs, self._df_tissue_positions, left_index=True, right_index=True)

    def on_data_loaded(self):
        self.merge_data()

    def load(self):
        self.read_img()
        self.read_adata()
        self.read_positions()
        return self.img, self.adata, self.df_tissue_positions
    
    def get_available_resolutions_data(self):
        results = {}
        if self.dir_base is None:
            return results

        for bs in self.get_available_sizes():
            adata_path, positions_path = self.get_default_path(bs, data='both')
            dataset = DataSet(
                img_path=self.img_path,
                adata_path=adata_path,
                positions_path=positions_path,
                bins_size=bs
            )
            dataset.read_adata()
            dataset.read_positions()
            results[bs] = (
                dataset.adata,
                dataset.df_tissue_positions
            )

        return results
