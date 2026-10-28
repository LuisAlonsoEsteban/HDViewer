
from ._sample import Sample

import math

class Spots(Sample):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_unit_area_um2(self):
        return math.pi * (55.0 / 2)**2