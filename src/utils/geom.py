
import numpy as np
import shapely
from shapely.geometry import Polygon


def exterior_coords(geoms):
    """Exterior ring of each polygon as two lists of vertex lists (xs, ys),
    vectorized; a geometry without exterior gets empty lists."""
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


def create_square(point, size=1):
    """Create a square polygon around the point with the given size.
    
    Args:
    point (shapely.geometry.Point): The central point of the square.
    size (float): The length of each side of the square.
    
    Returns:
    shapely.geometry.Polygon: The square polygon.
    """
    x, y = point.x, point.y
    half_size = size / 2
    square = Polygon([
        (x - half_size, y - half_size),
        (x - half_size, y + half_size),
        (x + half_size, y + half_size),
        (x + half_size, y - half_size)
    ])
    return square