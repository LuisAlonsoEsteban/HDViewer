
import os

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.load_data import DataSet
from src.tiles import ImageTiles

DATA_PATH = "/home/lalonso/PhD/VisHD/data/uploaded"

class DataInputPaths:

    def __init__(self, on_data_loaded_callback, on_image_loaded_callback):
        self._adata = None
        self._positions = None
        self._bins_size = None
        self._img = None

        self._id = None

        self.on_data_loaded_callback = on_data_loaded_callback
        self.on_image_loaded_callback = on_image_loaded_callback

    def write_file(self, path, data):
        with open(path, 'wb') as file:
            file.write(data.read())

    @property
    def id(self):
        if self._id:
            return self._id
        else:
            i = 1
            while True:
                candidate = str(i)
                if candidate not in os.listdir(DATA_PATH):
                    self._id = candidate
                    os.mkdir(os.path.join(DATA_PATH, self._id))
                    return self._id
                i += 1

    def is_data_loaded(self):
        return self._adata is not None and self._positions is not None and self._bins_size is not None
    
    def on_data_loaded(self):
        dataset = DataSet(
                    adata_path=os.path.join(DATA_PATH, self.id, "adata.h5ad"), 
                    positions_path=os.path.join(DATA_PATH, self.id, "positions.parquet"),
                    img_path="/",
                    bins_size=self.bins_size
                )
        dataset.read_adata()
        dataset.read_positions()

        self.on_data_loaded(dataset.adata, dataset.df_tissue_positions)

        return dataset.adata, dataset.df_tissue_positions

    def on_img_loaded(self):
        dataset = DataSet(
                    adata_path="/",
                    positions_path="/",
                    img_path=os.path.join(DATA_PATH, self.id, "image.btf"),
                    bins_size=2
                )
        dataset.read_img()

        self.on_image_loaded_callback(
            ImageTiles(dataset.img, os.path.join(self.img, "./tiles/"), segm_path=os.path.join(self.img, "./StarDist/"))
        )

        return dataset.img

    @property
    def adata(self):
        return self._adata

    @adata.setter
    def adata(self, value):
        self._adata = value
        self.write_file(os.path.join(DATA_PATH, self.id, "adata.h5ad"), value)
        if self.is_data_loaded():
            self.on_data_loaded()
            
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, value):
        self._positions = value
        self.write_file(os.path.join(DATA_PATH, self.id, "positions.parquet"), value)
        if self.is_data_loaded():
            self.on_data_loaded()

    @property
    def bins_size(self):
        return self._bins_size

    @bins_size.setter
    def bins_size(self, value):
        self._bins_size = value
        if self.is_data_loaded():
            self.on_data_loaded()

    @property
    def img(self):
        return self._img
    
    @img.setter
    def img(self, value):
        self._img = value
        self.write_file(os.path.join(DATA_PATH, self.id, "image.btf"), value)
        self.on_img_loaded()

