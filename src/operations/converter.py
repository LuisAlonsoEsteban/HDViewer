

import numpy as np
import pandas as pd

import math

class Converter:

    class VisiumHD:

        def __init__(self, pixels_per_bin=None, df_tissue_positions=None, bins_size=2.0, *, grid_fit_max_points=4096):
            self.grid_fit_max_points = grid_fit_max_points
            self.bins_size = bins_size
            self.pixels_per_bin = pixels_per_bin or self.get_pixels_per_bin(df_tissue_positions) if df_tissue_positions is not None else None
        
        def get_pixels_per_bin(self, df_tissue_positions:pd.DataFrame, max_iters=None):
            """Full-res pixels between neighbouring bins, from a least-squares fit of
            pixel positions against array positions (robust to rotation and flips)."""
            positions = df_tissue_positions
            sampled = positions.iloc[np.linspace(0, len(positions) - 1, min(self.grid_fit_max_points, len(positions)), dtype=int)]
            design = np.column_stack((sampled["array_col"], sampled["array_row"], np.ones(len(sampled))))
            pixels = sampled[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
            basis, _, rank, _ = np.linalg.lstsq(design, pixels, rcond=None)
            if rank < 3:
                return None
            return float(np.linalg.norm(basis[:2], axis=1).mean())

        def get_rotation(self, df_tissue_positions:pd.DataFrame, N=10, rad=True):
            rotations = []
            for i in range(N):
                fixed_dir = 'row' if i % 2 == 0 else 'col'
                orth_dir = 'row' if fixed_dir == 'col' else 'col'

                rand_index = np.random.randint(df_tissue_positions[f"array_{fixed_dir}"].min(), df_tissue_positions[f"array_{fixed_dir}"].max())
                sorted_rand_ax = df_tissue_positions[df_tissue_positions[f"array_{fixed_dir}"] == rand_index].sort_values(by=f"array_{orth_dir}")
                
                smaller_cat = (sorted_rand_ax.iloc[-1,:][f"pxl_{fixed_dir}_in_fullres"] - sorted_rand_ax.iloc[0,:][f"pxl_{fixed_dir}_in_fullres"]) * (1 if fixed_dir == 'col' else -1) # Vertical setup returns opposite angle, due to the sign switch of the smaller catete (a positive horizontal small catet turns into negative under vertical construction). This is to correct that.
                greater_cat = sorted_rand_ax.iloc[-1,:][f"pxl_{orth_dir}_in_fullres"] - sorted_rand_ax.iloc[0,:][f"pxl_{orth_dir}_in_fullres"]

                # Note that the direction of 'array_row' and 'pxl_row_in_fullres' are switched: as one increases, the other decreases and viceversa
                if fixed_dir == 'row':
                    smaller_cat *= -1
                else:
                    greater_cat *= -1
                    
                rot = math.atan2(smaller_cat, greater_cat)

                if abs(rot) > math.pi / 4:
                    rot = rot + math.pi / 2 * (1 if rot < 0 else -1)

                if not rad:
                    rot = math.degrees(rot)
                
                rotations.append(rot)

            return sum(rotations) / len(rotations)
        
        def get_length(self, n_bins:int):
            return n_bins * self.bins_size
        
        def get_bins(self, length:float):
            return length / self.bins_size
        
        def get_pixels(self, length:float):
            return self.get_bins(length) * self.pixels_per_bin

    class VisiumV2:
        
        SPOT_DIAMETER = 55.0
        CENTRE2CENTRE = 100.0

        def __init__(self):
            self.col_distance = self.CENTRE2CENTRE / 2
            self.row_distance = self.get_row_distance()

        def get_row_distance(self):
            return math.sqrt((self.CENTRE2CENTRE ** 2) - ((self.CENTRE2CENTRE / 2) ** 2))
        
        def get_length(self, n_spots:int, axis=0):
            return n_spots * (self.row_distance if axis == 0 else self.col_distance)
        
        def get_spots(self, length:float, axis=0):
            return length / (self.row_distance if axis == 0 else self.col_distance)
        

    def __init__(self, pixels_per_bin=None, df_tissue_positions:pd.DataFrame=None, bins_size=2.0, *, grid_fit_max_points=4096):
        self.visiumHD = self.VisiumHD(pixels_per_bin=pixels_per_bin, df_tissue_positions=df_tissue_positions, bins_size=bins_size, grid_fit_max_points=grid_fit_max_points)
        self.visiumV2 = self.VisiumV2()

    def bins2spots(self, n_bins:int, axis=0, shift=0):
        return self.visiumV2.get_spots(self.visiumHD.get_length(n_bins) + shift, axis)
    
    def spots2bins(self, spots:float, axis=0, shift=0):
        return self.visiumHD.get_bins(self.visiumV2.get_length(spots, axis) + shift)
    
    def spots2pixels(self, spots:float, axis=0, shift=0):
        return self.visiumHD.get_pixels(self.visiumV2.get_length(spots, axis) + shift)

    def get_bins(self, length:float):
        return self.visiumHD.get_bins(length)
    
    def get_spots(self, length:float, axis=0):
        return self.visiumV2.get_spots(length, axis)
    
    def get_pixels(self, length:float):
        return self.visiumHD.get_pixels(length)
