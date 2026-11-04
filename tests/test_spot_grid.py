import math
import unittest

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from hd_viewer.startup.operations.bins2spots import Bins2Spots
from hd_viewer.startup.operations.converter import Converter

PITCH = 58.4                      # full-res pixels per 16 µm bin
BIN_UM = 16


def positions(n=60, origin=(-3000.0, -2000.0), angle_deg=2.0):
    """A capture area starting at negative pixels, rotated, with rows running upwards
    in the array but downwards in the image (as in Space Ranger)."""
    a = math.radians(angle_deg)
    step_col = PITCH * np.array([math.cos(a), math.sin(a)])
    step_row = PITCH * np.array([math.sin(a), -math.cos(a)])      # flipped rows
    rows, cols = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    rows, cols = rows.ravel(), cols.ravel()
    xy = np.array(origin) + cols[:, None] * step_col + rows[:, None] * step_row
    names = [f"s_{r}_{c}" for r, c in zip(rows, cols)]
    df = pd.DataFrame({"in_tissue": 1, "array_row": rows, "array_col": cols,
                       "pxl_row_in_fullres": xy[:, 1], "pxl_col_in_fullres": xy[:, 0]}, index=names)
    df["index"] = df.index
    return df


def make_b2s(scale=(1.0, 1.34)):
    pos = positions()
    adata = ad.AnnData(np.ones((len(pos), 3)), obs=pd.DataFrame(index=pos.index),
                       var=pd.DataFrame(index=["a", "b", "c"]))
    return Bins2Spots(adata, pos, bins_size=BIN_UM, pixels_per_bin=None,
                      tile_scale_xfactor=scale[0], tile_scale_yfactor=scale[1])


class SpotGridTests(unittest.TestCase):
    def test_bin_pitch_estimate_survives_rotation_and_flips(self):
        self.assertAlmostEqual(Converter(df_tissue_positions=positions(), bins_size=16)
                               .visiumHD.pixels_per_bin, PITCH, places=6)

    def test_spots_cover_the_bins_not_the_image_origin(self):
        b2s = make_b2s()
        bins, spots = b2s.bin_squares.total_bounds, b2s.spot_polygons.total_bounds
        # Within one spot diameter of the bin grid on every side (was anchored at 0, 0)
        tolerance = 55 / BIN_UM * PITCH * 1.34
        for got, expected in zip(spots, bins):
            self.assertLess(abs(got - expected), tolerance)

    def test_buffer_radius_in_aligned_plot_units(self):
        import shapely
        from hd_viewer.startup.data.image.transform import SpatialTransform

        b2s = make_b2s()
        transform = SpatialTransform.fit(b2s.df_tissue_positions, (4000, 5000, 3), image_scale=.5)
        b2s.spatial_transform = transform
        spots = b2s.spot_polygons
        pixel_centres = spots[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy()
        centres = transform.points(pixel_centres)
        expected_radius = 27.5 / BIN_UM * transform.pitch
        coordinates, owners = shapely.get_coordinates(spots.geometry.to_numpy(), return_index=True)
        distances = np.linalg.norm(coordinates - centres[owners], axis=1)
        np.testing.assert_allclose(distances, expected_radius, atol=1e-9)
        np.testing.assert_allclose(spots.geometry.centroid.x, centres[:, 0], atol=1e-9)
        np.testing.assert_allclose(spots.geometry.centroid.y, centres[:, 1], atol=1e-9)

    def test_buffering_empty_spot_grid(self):
        b2s = Bins2Spots(None, positions(n=2), bins_size=BIN_UM)
        self.assertTrue(b2s.spot_polygons.empty)

    def test_visium_geometry_and_bins_inside_each_spot(self):
        b2s = make_b2s(scale=(1.0, 1.0))
        spots = b2s.spot_polygons
        centres = np.column_stack([spots.pxl_col_in_fullres, spots.pxl_row_in_fullres])
        nearest, _ = cKDTree(centres).query(centres, k=2)
        self.assertAlmostEqual(np.median(nearest[:, 1]) / PITCH * BIN_UM, 100.0, places=3)
        diameter = 2 * math.sqrt(spots.geometry.area.median() / math.pi) / PITCH * BIN_UM
        self.assertAlmostEqual(diameter, 55, delta=0.5)
        inside = gpd.sjoin(b2s.bin_squares[["geometry"]], spots[["geometry"]], predicate="within")
        per_spot = inside.groupby(inside.columns.drop("geometry")[0]).size()
        self.assertEqual(len(per_spot), len(spots))            # every spot is over bins
        self.assertTrue(per_spot.between(2, 6).all())


if __name__ == "__main__":
    unittest.main()
