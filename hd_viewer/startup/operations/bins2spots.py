

import pandas as pd
import geopandas as gpd
import numpy as np
import shapely

import scanpy as sc

from .converter import Converter
from .pooler import Pooler

from ..data.image.transform import fit_grid_basis


class Bins2Spots:

    def __init__(self, 
                 adata, 
                 tissue_positions, 
                 bins_size=None, 
                 pixels_per_bin:float=None,
                 tile_scale_xfactor=1.0, 
                 tile_scale_yfactor=1.0,
                 progress_callback=None,
                 spatial_transform=None,
                 *, grid_fit_max_points=4096, grid_max_condition_number=1e10,
                 spot_buffer_quad_segs=16, normalization_target_sum=1e4,
                 ):
        
        self.grid_fit_max_points = grid_fit_max_points
        self.grid_max_condition_number = grid_max_condition_number
        self.spot_buffer_quad_segs = spot_buffer_quad_segs
        self.normalization_target_sum = normalization_target_sum
        self.spatial_transform = spatial_transform
        self.adata_hd = adata
        self.df_tissue_positions = tissue_positions

        self.bins_size = bins_size
        self.converter = Converter(pixels_per_bin=pixels_per_bin, df_tissue_positions=tissue_positions, bins_size=bins_size, grid_fit_max_points=grid_fit_max_points)

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
    
    def grid_basis(self):
        if getattr(self, "_grid_basis", None) is not None:
            return self._grid_basis

        if self.spatial_transform is not None:
            self._grid_basis = (self.spatial_transform.inverse @ self.spatial_transform.grid_forward)[:2].T
            return self._grid_basis

        self._grid_basis = fit_grid_basis(self.df_tissue_positions, max_points=self.grid_fit_max_points,
                                          max_condition_number=self.grid_max_condition_number)
        return self._grid_basis

    @property
    def display_scale(self):
        """Full-res pixels -> figure units, per axis (the image is shown as a square)."""
        return np.array([self.tile_scale_xfactor, self.tile_scale_yfactor])

    def create_bin_squares(self, subset=False):
        positions = self.df_tissue_positions
        centers = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        if self.spatial_transform is not None:
            centers = self.spatial_transform.points(centers)
            geometry = self.spatial_transform.rectangles(centers)
        else:
            basis = self.grid_basis()
            offsets = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ basis[:2]
            corners = (centers[:, None, :] + offsets) * self.display_scale
            geometry = shapely.polygons(corners)
        self._bin_squares = gpd.GeoDataFrame(self.df_tissue_positions, geometry=geometry)
        return self._bin_squares[self._bin_squares["in_tissue"] == 1] if subset else self._bin_squares

    def create_spot_polygons(self):
        """Visium v2-like spots (55 µm, 100 µm centre to centre, hexagonal) over the bins.

        The grid is laid out in the bins' array coordinates (in bins) and mapped to
        pixels with the same fitted transform as the bins, so it shares their origin,
        rotation and orientation, then scaled like them for display.
        """

        v2 = self.converter.visiumV2
        bin_um = float(self.bins_size)
        positions = self.df_tissue_positions
        col_min, row_min = positions["array_col"].min(), positions["array_row"].min()
        nbins_cols = positions["array_col"].max() - col_min + 1
        nbins_rows = positions["array_row"].max() - row_min + 1

        col_step, row_step = v2.col_distance / bin_um, v2.row_distance / bin_um
        radius = v2.SPOT_DIAMETER / 2 / bin_um
        nspots_cols = int(nbins_cols * bin_um / v2.col_distance)
        nspots_rows = int(nbins_rows * bin_um / v2.row_distance)

        grid_rows, grid_cols = [], []
        for y in range(nspots_rows):
            for x in range(y % 2, nspots_cols, 2):
                grid_rows.append(y)
                grid_cols.append(x)
            if self.progress_callback is not None:
                self.progress_callback(1, y + 1, nspots_rows)
        grid_rows, grid_cols = np.array(grid_rows), np.array(grid_cols)

        # Centres in array coordinates (bin units; bin 0 spans -0.5..0.5)
        centre_col = col_min - 0.5 + radius + grid_cols * col_step
        centre_row = row_min - 0.5 + radius + grid_rows * row_step

        basis = self.grid_basis()
        def to_pixels(cols, rows):
            design = np.stack([cols, rows, np.ones_like(cols)], axis=-1)
            return design @ basis

        centres = to_pixels(centre_col, centre_row)

        if self.spatial_transform is not None:
            display_centres = self.spatial_transform.points(centres)
            geometries = shapely.buffer(shapely.points(display_centres),
                                        radius * self.spatial_transform.pitch, quad_segs=self.spot_buffer_quad_segs)
        else:
            seeds = shapely.points(np.column_stack((centre_col, centre_row)))
            geometries = shapely.buffer(seeds, radius, quad_segs=self.spot_buffer_quad_segs)
            geometries = shapely.transform(geometries, lambda xy: (xy @ basis[:2] + basis[2]) * self.display_scale)

        metadata = {
            'index': np.arange(len(grid_rows)),
            'array_row': grid_rows,
            'array_col': grid_cols,
            # Full-resolution pixel position of each spot centre (as in Space Ranger)
            'pxl_row_in_fullres': centres[:, 1],
            'pxl_col_in_fullres': centres[:, 0],
        }

        self._spot_polygons = gpd.GeoDataFrame(metadata, geometry=geometries)
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
            sc.pp.normalize_total(self.adata_hd, target_sum=self.normalization_target_sum)
            sc.pp.normalize_total(self.adata_v2, target_sum=self.normalization_target_sum)

        return self.adata_v2
