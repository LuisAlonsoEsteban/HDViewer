import re

from fastapi import HTTPException, Request, Response

from ... import settings


def current_session(request: Request, response: Response):
    """FastAPI dependency: the caller's session (created on first use). Only the
    entry points of a visit use it (GET/POST /api/load, uploads)."""
    return request.app.state.sessions.resolve(request, response)


def existing_session(request: Request):
    """FastAPI dependency: the caller's session; never creates one, so requests
    without a valid cookie (an expired session, a crawler) do not use a slot."""
    session = request.app.state.sessions.existing(request)
    if session is None:
        raise HTTPException(409, "No dataset is loaded.")
    return session


def public_message(exc):
    """An exception's message for the browser, without server folders (session
    folders, the data and package directories)."""
    message = str(exc)
    for folder in (settings.UPLOAD_DIR, settings.DATA_DIR, settings.BASE_DIR.parent):
        message = re.sub(re.escape(str(folder)) + r"/(?:[0-9a-f]{32}/)?", "", message)
    return message
