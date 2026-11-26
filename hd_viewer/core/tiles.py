"""Tile lookup for the image tiles route of the HTTP API.

Tile URLs are /<key>/<zoom>/<x>/<y>[.<ext>]; keys are content-versioned
(see startup.data.image.image.register_tiles), so a found tile never changes.
"""
import os

from .. import settings
from ..startup.data.image.image import (
    TILE_CACHE_CONTROL, TILE_MIMETYPES, TILE_NAME, _TILE_LOCK, _TILE_SETS,
)

__all__ = ["resolve_tile", "tile_cache_control", "forget_tiles_under", "TILE_CACHE_CONTROL"]

# Tiles of uploaded images belong to one visitor: shared caches must not keep them
PRIVATE_TILE_CACHE_CONTROL = 'private, max-age=31536000, immutable'


def _under(path, folder):
    folder = os.path.join(os.path.abspath(str(folder)), "")
    return os.path.abspath(str(path)).startswith(folder)


def tile_cache_control(path):
    return PRIVATE_TILE_CACHE_CONTROL if _under(path, settings.UPLOAD_DIR) else TILE_CACHE_CONTROL


def forget_tiles_under(folder):
    """Stop serving the tile sets stored in `folder` (a closed session's folder)."""
    with _TILE_LOCK:
        for key in [k for k, (p, _) in _TILE_SETS.items() if _under(p, folder)]:
            del _TILE_SETS[key]


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
