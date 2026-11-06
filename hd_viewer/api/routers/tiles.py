from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, Response

from ...core.tiles import TILE_CACHE_CONTROL, resolve_tile

router = APIRouter()
CORS = {"Access-Control-Allow-Origin": "*"}


@router.get("/tiles/{key}/{zoom}/{x}/{tail}")
def tile(key: str, zoom: int, x: int, tail: str, request: Request):
    """Image tiles; URLs are content-versioned, so found tiles are cached forever."""
    y, _, ext = tail.partition(".")
    try:
        found = resolve_tile(key, zoom, x, int(y), ext or None)
    except ValueError:
        found = None
    if found is None:
        raise HTTPException(404, "Tile not found", headers={**CORS, "Cache-Control": "no-store"})
    path, mimetype, etag = found
    etag = f'"{etag}"'
    headers = {**CORS, "Cache-Control": TILE_CACHE_CONTROL, "ETag": etag}
    if request.headers.get("if-none-match") == etag:
        return Response(status_code=304, headers=headers)
    return FileResponse(path, media_type=mimetype, headers=headers)
