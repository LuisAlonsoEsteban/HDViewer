
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely

import scanpy as sc

from .converter import Converter
from .pooler import Pooler


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

    def grid_basis(self):
        """Affine map (array_col, array_row, 1) -> full-res pixels (col, row), fitted to
        the positions. Captures the grid's origin, pitch, rotation and the flipped row
        direction; bins and pseudospots both go through it, so they always line up."""
        if getattr(self, "_grid_basis", None) is not None:
            return self._grid_basis
        positions = self.df_tissue_positions
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
            origin = pixel_centers.mean(axis=0) if len(pixel_centers) else np.zeros(2)
            basis = np.array([[pitch, 0], [0, pitch], origin])
        self._grid_basis = basis
        return basis

    @property
    def display_scale(self):
        """Full-res pixels -> figure units, per axis (the image is shown as a square)."""
        return np.array([self.tile_scale_xfactor, self.tile_scale_yfactor])

    def create_bin_squares(self, subset=False):
        positions = self.df_tissue_positions
        # A reflected or rotated array cannot be represented by swapping rectangle
        # dimensions: offset the corners along the fitted grid directions.
        basis = self.grid_basis()
        offsets = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ basis[:2]
        centers = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        corners = (centers[:, None, :] + offsets) * self.display_scale
        geometry = shapely.polygons(corners)
        self._bin_squares = gpd.GeoDataFrame(self.df_tissue_positions, geometry=geometry)
        return self._bin_squares[self._bin_squares["in_tissue"] == 1] if subset else self._bin_squares

    # Vertices of each pseudospot outline
    SPOT_VERTICES = 32

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

        # Spacing in bins: neighbouring columns 50 µm apart (alternate rows shifted),
        # rows sqrt(100² - 50²) µm apart; the first spot sits one radius from the edge
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

        angles = np.linspace(0, 2 * np.pi, self.SPOT_VERTICES, endpoint=False)
        outline_cols = centre_col[:, None] + radius * np.cos(angles)
        outline_rows = centre_row[:, None] + radius * np.sin(angles)
        outlines = to_pixels(outline_cols, outline_rows) * self.display_scale
        centres = to_pixels(centre_col, centre_row)

        metadata = {
            'index': np.arange(len(grid_rows)),
            'array_row': grid_rows,
            'array_col': grid_cols,
            # Full-resolution pixel position of each spot centre (as in Space Ranger)
            'pxl_row_in_fullres': centres[:, 1],
            'pxl_col_in_fullres': centres[:, 0],
        }
        self._spot_polygons = gpd.GeoDataFrame(metadata, geometry=shapely.polygons(outlines))
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
