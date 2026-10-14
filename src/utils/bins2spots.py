
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely

from shapely.geometry import Point

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
                 tile_scale_yfactor=1.0,
                 progress_callback=None,
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

        self.progress_callback = progress_callback


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
        positions = self.df_tissue_positions
        # Fit the two grid directions in image coordinates. A reflected or
        # rotated array cannot be represented by swapping rectangle dimensions.
        # Sampling bounds the cost for high-resolution position tables.
        sampled = positions.iloc[np.linspace(
            0, len(positions) - 1, min(4096, len(positions)), dtype=int
        )]
        design = np.column_stack((sampled["array_col"], sampled["array_row"],
                                  np.ones(len(sampled))))
        pixel_centers = sampled[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        basis, _, rank, _ = np.linalg.lstsq(design, pixel_centers, rcond=None)
        if rank < 3:
            # A single row/column cannot determine both directions.
            pitch = self.converter.visiumHD.pixels_per_bin
            basis = np.array([[pitch, 0], [0, pitch]])
        offsets = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ basis[:2]
        centers = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        # Rotate in raw image coordinates, then scale each display axis.
        corners = (centers[:, None, :] + offsets) * np.array(
            [self.tile_scale_xfactor, self.tile_scale_yfactor]
        )
        geometry = shapely.polygons(corners)
        self._bin_squares = gpd.GeoDataFrame(self.df_tissue_positions, geometry=geometry)
        return self._bin_squares[self._bin_squares["in_tissue"] == 1] if subset else self._bin_squares

    def create_spot_polygons(self):
        spots = []
        metadata = {'index': [], 'array_row': [], 'array_col': [], 'pxl_row_in_fullres': [], 'pxl_col_in_fullres': []}
        
        nbins_rows, nbins_cols = pd.unique(self.bin_squares["array_row"]).shape[0], pd.unique(self.bin_squares["array_col"]).shape[0] 
        nspots_rows = int(self.converter.bins2spots(nbins_rows, axis=0))
        nspots_cols = int(self.converter.bins2spots(nbins_cols, axis=1))

        id = 0
        scaled_radius = self.converter.get_pixels(self.converter.visiumV2.SPOT_DIAMETER / 2)
        max_scaled_y = self.converter.spots2pixels(nspots_rows)

        for i, y in tqdm(enumerate(range(nspots_rows))):
            scaled_y = self.converter.spots2pixels(y, axis=0, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)
            for x in range(i % 2, nspots_cols, 2):
                metadata['index'].append(id)
                id += 1

                metadata['array_row'].append(y)
                metadata['array_col'].append(x)

                scaled_x = self.converter.spots2pixels(x, axis=1, shift=self.converter.visiumV2.SPOT_DIAMETER / 2)

                # metadata['pxl_row_in_fullres'].append(scaled_x)
                metadata['pxl_row_in_fullres'].append(self.tile_scale_yfactor * (max_scaled_y - scaled_y)) # # Note that the direction of 'array_row' and 'pxl_row_in_fullres' are switched: as one increases, the other decreases and viceversa
                metadata['pxl_col_in_fullres'].append(self.tile_scale_xfactor * scaled_x)

                spots.append(Point(scaled_x, scaled_y).buffer(scaled_radius))

            if self.progress_callback is not None:
                self.progress_callback(1, i + 1, nspots_rows)

        self._spot_polygons = gpd.GeoDataFrame(metadata, geometry=spots)
        self._spot_polygons.index.name = 'spot_id'  # avoid conflict with 'index' column in reset_index()

        return self._spot_polygons
    
    def convert(self, normalize=True):
        pooler = Pooler(
            self.bin_squares,
            self.adata_hd,
            self.spot_polygons,
            keep_cols=['array_row', 'array_col', 'pxl_row_in_fullres','pxl_col_in_fullres']
        )

        self.adata_v2 = pooler.run(predicate='within')

        if normalize:
            sc.pp.normalize_total(self.adata_hd, target_sum=1e4)
            sc.pp.normalize_total(self.adata_v2, target_sum=1e4)

        return self.adata_v2
