

from ._sample import Sample
from ..utils.geom import exterior_coords

import pandas as pd
import numpy as np
import shapely
import math

class Cells(Sample):
    # Levels of detail: level k draws the polygons simplified with tolerance
    # LOD_BASE_TOLERANCE * 2**k (figure units); None draws the original polygons
    LOD_PREFIX = 'lod:'
    LOD_BASE_TOLERANCE = 1.0
    LOD_MAX_LEVEL = 8
    LOD_TOLERANCE_PX = 1.0
    LOD_HYSTERESIS = 0.15

    def __init__(self, *args, lod_prefix=None, lod_base_tolerance=None, lod_max_level=None,
                 lod_tolerance_px=None, lod_hysteresis=None, **kwargs):
        super().__init__(*args, **kwargs)
        for name, value in (('LOD_PREFIX', lod_prefix), ('LOD_BASE_TOLERANCE', lod_base_tolerance),
                            ('LOD_MAX_LEVEL', lod_max_level), ('LOD_TOLERANCE_PX', lod_tolerance_px),
                            ('LOD_HYSTERESIS', lod_hysteresis)):
            if value is not None:
                setattr(self, name, value)
        self._simplified = {}

    def get_lod_level(self, pixel_size, current=None):
        tolerance = pixel_size * self.LOD_TOLERANCE_PX
        if not math.isfinite(tolerance) or tolerance <= 0:
            return current

        value = math.log2(tolerance / self.LOD_BASE_TOLERANCE)
        k = self.get_lod_index(current)
        if not (k - self.LOD_HYSTERESIS <= value < k + 1 + self.LOD_HYSTERESIS):
            k = max(-1, min(math.floor(value), self.LOD_MAX_LEVEL))

        return None if k < 0 else f"{self.LOD_PREFIX}{k}"

    def get_lod_index(self, resolution_level):
        if isinstance(resolution_level, str) and resolution_level.startswith(self.LOD_PREFIX):
            return int(resolution_level[len(self.LOD_PREFIX):])
        return -1

    def get_tolerance(self, resolution_level):
        k = self.get_lod_index(resolution_level)
        return 0.0 if k < 0 else self.LOD_BASE_TOLERANCE * 2 ** k

    def get_simplified(self, tolerance):
        if tolerance not in self._simplified:
            geoms = np.asarray(self.gdf.geometry.values, dtype=object)
            self._simplified[tolerance] = shapely.simplify(geoms, tolerance, preserve_topology=True)
        return self._simplified[tolerance]

    def get_coordinates(self, gdf, tolerance=0.0):
        if tolerance > 0:
            geoms = self.get_simplified(tolerance)[self.gdf.index.get_indexer(gdf.index)]
        else:
            geoms = np.asarray(gdf.geometry.values, dtype=object)
            
        xs, ys = exterior_coords(geoms)
        return pd.DataFrame({"x": xs, "y": ys}, index=gdf.index)

    def fetch_data(self, bbox=None, colnames=None, return_coords=True, resolution_level=None) -> pd.DataFrame:
        return super().fetch_data(bbox=bbox, colnames=colnames, return_coords=return_coords,
                                  tolerance=self.get_tolerance(resolution_level))