
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point

import scipy.sparse as sparse

import scanpy as sc

from tqdm import tqdm

from .converter import Converter
from .geom import create_square
from .operations import Pooler


class Bins2Spots:

    def __init__(self, 
                 adata, 
                 tissue_positions, 
                 bins_size=None, 
                 pixels_per_bin:float=None,
                 tile_scale_xfactor=1.0, 
                 tile_scale_yfactor=1.0
                 ):
        
        self.adata_hd = adata
        self.df_tissue_positions = tissue_positions

        self.bins_size = bins_size
        self.converter = Converter(pixels_per_bin=pixels_per_bin, df_tissue_positions=tissue_positions, bins_size=bins_size)

        self._bin_squares = None
        self._spot_polygons = None
        self.adata_v2 = None

        self.tile_scale_xfactor = tile_scale_xfactor
        self.tile_scale_yfactor = tile_scale_yfactor


    @property
    def bin_squares(self):
        if self._bin_squares is None:
            self.create_bin_squares()
        return self._bin_squares
    
    @property
    def spot_polygons(self):
        if self._spot_polygons is None:
            self.create_spot_polygons()
        return self._spot_polygons
    
    def get_rotation(self, rad=False):
        return self.converter.visiumHD.get_rotation(self.df_tissue_positions, rad=rad)

    def create_bin_squares(self, subset=False):
        geometry = [create_square(Point(y * self.tile_scale_xfactor, x * self.tile_scale_yfactor), self.converter.visiumHD.pixels_per_bin) 
                    for (x,y) in tqdm(zip(self.df_tissue_positions['pxl_row_in_fullres'], self.df_tissue_positions['pxl_col_in_fullres']))]
        self._bin_squares = gpd.GeoDataFrame(self.df_tissue_positions, geometry=geometry)
        return self._bin_squares[self._bin_squares["in_tissue"]] if subset else self._bin_squares

    def create_spot_polygons(self):
        spots = []
        metadata = {'index': [], 'array_row': [], 'array_col': [], 'pxl_row_in_fullres': [], 'pxl_col_in_fullres': []}
        
        nbins_rows, nbins_cols = pd.unique(self.bin_squares["array_row"]).shape[0], pd.unique(self.bin_squares["array_col"]).shape[0] 
        nspots_rows = int(self.converter.bins2spots(nbins_rows, axis=0))
        nspots_cols = int(self.converter.bins2spots(nbins_cols, axis=1))

        id = 0
        scaled_radius = self.converter.get_pixels(self.converter.visiumV2.SPOT_DIAMETER / 2)
        max_scaled_x = self.converter.spots2pixels(nspots_rows)
        for i, x in enumerate(tqdm(range(nspots_rows))):
            for y in range(i % 2, nspots_cols, 2):
                metadata['index'].append(id)
                id += 1

                metadata['array_row'].append(x)
                metadata['array_col'].append(y)

                scaled_x = self.converter.spots2pixels(x, axis=0, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)
                scaled_y = self.converter.spots2pixels(y, axis=1, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)
                # metadata['pxl_row_in_fullres'].append(scaled_x)
                metadata['pxl_row_in_fullres'].append(self.tile_scale_xfactor * (max_scaled_x - scaled_x)) # # Note that the direction of 'array_row' and 'pxl_row_in_fullres' are switched: as one increases, the other decreases and viceversa
                metadata['pxl_col_in_fullres'].append(self.tile_scale_yfactor * scaled_y)

                spots.append(Point(scaled_y, scaled_x).buffer(scaled_radius))

        self._spot_polygons = gpd.GeoDataFrame(metadata, geometry=spots)

        return self._spot_polygons
    
    def convert(self, normalize=True):
        pooler = Pooler(
            self.bin_squares, 
            self.adata_hd, 
            self.spot_polygons,
            keep_cols=['array_row', 'array_col', 'pxl_row_in_fullres','pxl_col_in_fullres']
        )
        self.adata_v2 = pooler.run(predicate='intersects')

        if normalize:
            sc.pp.normalize_total(self.adata_hd, target_sum=1e4)
            sc.pp.normalize_total(self.adata_v2, target_sum=1e4)
        
        return self.adata_v2 
