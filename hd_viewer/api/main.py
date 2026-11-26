"""HTTP API and static files of the React front end.

Run with ``uvicorn hd_viewer.api.main:app`` (one process: sessions, loaded
datasets and the tile registry are held in memory). Figures are still built with
Bokeh on the server and drawn by BokehJS in the browser.
"""
import asyncio
import contextlib
import os
from pathlib import Path

import bokeh
from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from ..core.session import purge_orphan_folders
from .sessions import SessionRegistry, TooManySessions, set_session_cookie

PACKAGE = Path(__file__).resolve().parents[1]
DIST = Path(os.environ.get("VHD_FRONTEND_DIST", Path(__file__).resolve().parent / "dist"))
BOKEHJS = Path(bokeh.__file__).resolve().parent / "server" / "static" / "js"
SWEEP_SECONDS = 30
CAPACITY_RETRY_SECONDS = 60      # when every session is in use, the browser retries after this


def create_app(registry=None, sweep_seconds=SWEEP_SECONDS, warm_up=None):
    from .routers import (analysis, comparison, config, datasets, de, deconvolution, embedding, events, layers, load,
                          metadata, plugins, regions, summary, tiles, uploads)

    registry = SessionRegistry() if registry is None else registry

    @contextlib.asynccontextmanager
    async def lifespan(app):
        # Sessions live in memory: every session folder on disk at startup is an orphan
        await asyncio.to_thread(purge_orphan_folders, registry.active_folders())
        if warm_up if warm_up is not None else os.environ.get("VHD_WARM_UP", "1") != "0":
            from ..core.warmup import warm_up as start_warm_up
            start_warm_up()

        async def sweep():
            while True:
                await asyncio.sleep(sweep_seconds)
                registry.sweep()
                await asyncio.to_thread(purge_orphan_folders, registry.active_folders())
        task = asyncio.create_task(sweep())
        try:
            yield
        finally:
            task.cancel()

    app = FastAPI(title="HDViewer API", lifespan=lifespan,
                  root_path=os.environ.get("VHD_ROOT_PATH", ""))
    app.state.sessions = registry

    @app.middleware("http")
    async def session_cookie(request: Request, call_next):
        """Set the cookie of a session created by this request on every response
        (also raw and error responses, which skip FastAPI's dependency headers)."""
        return set_session_cookie(request, await call_next(request))

    @app.exception_handler(TooManySessions)
    async def at_capacity(request, exc):
        """Every analysis session is in use: an expected, temporary state (not an
        error). The browser shows a notice and retries after `retry_after` seconds."""
        return JSONResponse(
            {"detail": ("All analysis sessions on this server are in use right now. Sessions free up as "
                        "other visitors finish, so please try again in a couple of minutes."),
             "reason": "capacity", "retry_after": CAPACITY_RETRY_SECONDS},
            status_code=503, headers={"Retry-After": str(CAPACITY_RETRY_SECONDS)})

    for router in (config, datasets, plugins, events, tiles, uploads, load, layers, analysis, regions, de,
                   deconvolution, metadata, embedding, summary, comparison):
        app.include_router(router.router, prefix="/api")

    # BokehJS from the installed Python package: the browser always runs the
    # same Bokeh version that generated the figures
    app.mount("/static/bokehjs", StaticFiles(directory=BOKEHJS), name="bokehjs")
    app.mount("/static/assets", StaticFiles(directory=PACKAGE / "static"), name="assets")

    @app.get("/{path:path}", include_in_schema=False)
    async def frontend(path: str, request: Request):
        """The built React app (single-page: unknown paths serve index.html)."""
        if path.startswith("api/"):
            return JSONResponse({"detail": "Not Found"}, status_code=404)
        if not (DIST / "index.html").exists():
            return PlainTextResponse("The front end is not built: run `npm run build` in frontend/.",
                                     status_code=503)
        target = (DIST / path).resolve()
        if path and target.is_file() and DIST.resolve() in target.parents:
            # Vite's assets/ files carry a content hash in their name
            immutable = path.startswith("assets/")
            return FileResponse(target, headers={"Cache-Control": "public, max-age=31536000, immutable"}
                                if immutable else None)
        return FileResponse(DIST / "index.html")

    return app


app = create_app()
