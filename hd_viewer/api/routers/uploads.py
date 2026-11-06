"""Chunked uploads over HTTP: POST announces a file, PUT sends raw-byte chunks.

Validation (names, extensions, size limit, disk space) and assembly are the
shared core.uploads.UploadField, so both front ends accept the same files.
"""
import secrets

from fastapi import APIRouter, Depends, HTTPException, Request
from pydantic import BaseModel

from ...core.uploads import UploadField
from . import current_session

router = APIRouter()

FIELDS = {
    "expression": ((".h5",), True),
    "positions": ((".parquet",), True),
    "image": ((".btf", ".png", ".tif", ".tiff"), True),
    "reference": ((".h5ad", ".h5"), False),
    "reference_annotations": ((".csv", ".parquet"), False),
}


class UploadStart(BaseModel):
    field: str
    filename: str
    size: int


def _field(session, key):
    if key not in FIELDS:
        raise HTTPException(404, f"Unknown upload field: {key}")
    if key not in session.uploads:
        session.uploads[key] = UploadField(key, FIELDS[key][0], lambda: session.directory)
    return session.uploads[key]


def describe(session):
    return {key: dict(status=field.status, filename=field.filename, size=field.size,
                      received=field.received, error=field.error, required=FIELDS[key][1])
            for key, field in ((key, _field(session, key)) for key in FIELDS)}


def _unlocked(session):
    if session.load.get("state") == "running":
        raise HTTPException(409, "Uploads are locked while data is loading.")


@router.get("/uploads")
def uploads(session=Depends(current_session)):
    return describe(session)


@router.post("/uploads")
def start_upload(body: UploadStart, session=Depends(current_session)):
    _unlocked(session)
    with session.lock:
        field = _field(session, body.field)
        upload_id = secrets.token_urlsafe(12)
        error = field.begin(upload_id, body.filename, body.size)
    if error is not None:
        raise HTTPException(422, error)
    return dict(upload_id=upload_id, chunk_bytes=field.chunk_bytes)


@router.put("/uploads/{key}/{upload_id}/{index}")
async def upload_chunk(key: str, upload_id: str, index: int, request: Request,
                       session=Depends(current_session)):
    _unlocked(session)
    field = _field(session, key)
    chunk = await request.body()
    with session.lock:
        if not field.is_current(upload_id):
            raise HTTPException(409, "This upload is no longer active.")
        error = field.write(index, chunk)
    if error is not None:
        raise HTTPException(422, error)
    return dict(received=field.received, size=field.size, complete=field.complete)


@router.delete("/uploads/{key}")
def remove_upload(key: str, session=Depends(current_session)):
    _unlocked(session)
    with session.lock:
        _field(session, key).clear()
    return describe(session)
