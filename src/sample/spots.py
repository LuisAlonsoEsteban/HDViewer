
from ..operations.converter import Converter
from ._sample import Sample

import math

class Spots(Sample):
    def __init__(self, *args, radius=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius

    def get_unit_area_um2(self):
        return math.pi * (Converter.VisiumV2.SPOT_DIAMETER / 2)**2

    def get_coordinates(self, gdf):
        coordinates = super().get_coordinates(gdf)
        if self.radius is not None:
            coordinates['radius'] = self.radius
        else:
            bounds = gdf.geometry.bounds
            coordinates['radius'] = (bounds.maxx - bounds.minx) / 2
        return coordinates
