
from shapely.geometry import Polygon

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