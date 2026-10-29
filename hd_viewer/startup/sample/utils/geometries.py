
import numpy as np
import shapely

def exterior_coords(geoms):
    geoms = np.asarray(geoms, dtype=object)
    coords, owner = shapely.get_coordinates(shapely.get_exterior_ring(geoms), return_index=True)
    counts = np.bincount(owner, minlength=len(geoms))

    if len(counts) and (counts == counts[0]).all() and counts[0] > 0:
        # Same number of vertices everywhere (e.g. bins): one reshape, no split
        rings = coords.reshape(len(geoms), counts[0], 2)
        return rings[:, :, 0].tolist(), rings[:, :, 1].tolist()
        
    bounds = np.cumsum(counts)[:-1]

    return ([v.tolist() for v in np.split(coords[:, 0], bounds)],
            [v.tolist() for v in np.split(coords[:, 1], bounds)])