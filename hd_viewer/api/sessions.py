"""Browser sessions of the HTTP API.

A session is identified by a random cookie (`hdv_sid`); its server folder is a
hash of it (core.session), so a client cannot choose a path. Sessions end when
their events socket has been closed for `GRACE_SECONDS` (a page reload keeps
them) or after `ttl` seconds without any request; ending one cancels its jobs
and removes its folder.
"""
import asyncio
import os
import secrets
import threading
import time
from dataclasses import dataclass, field

from fastapi import Request, Response

from ..core.session import cleanup_session, session_upload_dir

COOKIE = "hdv_sid"
GRACE_SECONDS = 60
DEFAULT_TTL = int(os.environ.get("VHD_SESSION_TTL", 1800))
MAX_SESSIONS = int(os.environ.get("VHD_MAX_SESSIONS", 1))


class TooManySessions(Exception):
    pass


@dataclass
class Session:
    sid: str
    last_seen: float = field(default_factory=time.monotonic)
    lock: threading.RLock = field(default_factory=threading.RLock)
    subscribers: set = field(default_factory=set)
    loop: asyncio.AbstractEventLoop | None = None
    uploads: dict = field(default_factory=dict)       # field key -> core.uploads.UploadField
    load: dict = field(default_factory=dict)          # load job status: state, tasks, error
    loaded: object = None                             # core.pipeline.LoadedDataset
    views: dict = field(default_factory=dict)         # layer id -> api.layers.LayerView
    comparison: dict = field(default_factory=dict)    # sample comparison: job, result (api.routers.comparison)
    disconnected_at: float | None = None

    @property
    def directory(self):
        return session_upload_dir(self.sid)

    def touch(self):
        self.last_seen = time.monotonic()

    def publish(self, event):
        """Send `event` to every open events socket; safe from any thread."""
        loop = self.loop
        if loop is None or loop.is_closed():
            return
        for queue in list(self.subscribers):
            loop.call_soon_threadsafe(queue.put_nowait, event)


class SessionRegistry:
    def __init__(self, ttl=DEFAULT_TTL, max_sessions=MAX_SESSIONS):
        self.ttl = ttl
        self.max_sessions = max_sessions
        self._sessions = {}
        self._lock = threading.Lock()

    def __len__(self):
        return len(self._sessions)

    def get(self, sid):
        with self._lock:
            return self._sessions.get(sid)

    def create(self):
        with self._lock:
            if len(self._sessions) >= self.max_sessions:
                raise TooManySessions()
            sid = secrets.token_urlsafe(32)
            session = self._sessions[sid] = Session(sid)
            return session

    def resolve(self, request: Request, response: Response):
        """The request's session, creating it (and its cookie) if needed."""
        session = self.get(request.cookies.get(COOKIE, ""))
        if session is None:
            session = self.create()
            response.set_cookie(COOKIE, session.sid, httponly=True, samesite="lax",
                                secure=request.url.scheme == "https", path="/")
        session.touch()
        return session

    def close(self, sid):
        with self._lock:
            session = self._sessions.pop(sid, None)
        if session is not None:
            threading.Thread(target=cleanup_session, args=(sid,), daemon=True).start()
        return session

    def sweep(self, now=None):
        """End idle sessions and those whose socket closed more than GRACE_SECONDS ago."""
        now = time.monotonic() if now is None else now
        with self._lock:
            expired = [sid for sid, session in self._sessions.items()
                       if not session.subscribers and (
                           now - session.last_seen > self.ttl
                           or (session.disconnected_at is not None
                               and now - session.disconnected_at > GRACE_SECONDS))]
        for sid in expired:
            self.close(sid)
        return expired
