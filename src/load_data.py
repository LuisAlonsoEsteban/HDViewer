
import os
import re

from tifffile import imread
import pandas as pd
import scanpy as sc


class DataSet:

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
        
        
        self.img_path = img_path or os.path.join(dir_base, img_file)

        self.adata_path = adata_path or os.path.join(dir_base, f"binned_outputs/square_{bins_size:03}um/filtered_feature_bc_matrix.h5")
        self.positions_path = positions_path or os.path.join(dir_base, f"binned_outputs/square_{bins_size:03}um/spatial/tissue_positions.parquet")

        if bins_size:
            self.bins_size = bins_size
        else:
            assert self.get_bins_size(adata_path) == self.get_bins_size(positions_path)
            self.bins_size = self.get_bins_size(adata_path)

    @staticmethod
    def get_bins_size(path:str) -> int:
        try:
            bins_size_comp = path.split(os.sep)[-2]
            return int(re.search(r'\d{3}', bins_size_comp).group())
        except IndexError:
            raise ValueError(f"Path {path} has no parent directory on which bins size can be retrieved.")
        except AttributeError:
            raise ValueError(f"Path {path} has no bins size in it. You must provide it.")
        
    def read_img(self):
        self.img = imread(self.img_path)

    def read_positions(self):
        self._df_tissue_positions = pd.read_parquet(self.positions_path)
        self._df_tissue_positions = self._df_tissue_positions.set_index('barcode')
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


