"""Tile lookup shared by every server that serves image tiles (Flask, FastAPI).

Tile URLs are /<key>/<zoom>/<x>/<y>[.<ext>]; keys are content-versioned
(see startup.data.image.image.register_tiles), so a found tile never changes.
"""
import os

from ..startup.data.image.image import (
    TILE_CACHE_CONTROL, TILE_MIMETYPES, TILE_NAME, _TILE_LOCK, _TILE_SETS,
)

__all__ = ["resolve_tile", "TILE_CACHE_CONTROL"]


def resolve_tile(key, zoom, x, y, ext=None):
    """(path, mimetype, etag) of a registered tile, or None when there is none."""
    with _TILE_LOCK:
        tile_set = _TILE_SETS.get(key)
    if tile_set is None:
        return None
    tiles_path, tile_ext = tile_set
    if ext is not None and ext != tile_ext:
        return None
    path = os.path.join(tiles_path, TILE_NAME.format(zoom=zoom, x=x, y=y, ext=tile_ext))
    if not os.path.exists(path):
        return None
    return path, TILE_MIMETYPES.get(tile_ext, 'application/octet-stream'), f"{key}-{zoom}-{x}-{y}"
