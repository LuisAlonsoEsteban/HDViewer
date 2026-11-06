"""HTTP API and static files of the React front end.

Run with ``uvicorn hd_viewer.api.main:app`` (one process: sessions and tiles are
held in memory). The Bokeh-server UI (`bokeh serve hd_viewer`) keeps working
side by side; both use the same core.
"""
import asyncio
import contextlib
import os
from pathlib import Path

import bokeh
from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from .sessions import SessionRegistry, TooManySessions

PACKAGE = Path(__file__).resolve().parents[1]
DIST = Path(os.environ.get("VHD_FRONTEND_DIST", Path(__file__).resolve().parent / "dist"))
BOKEHJS = Path(bokeh.__file__).resolve().parent / "server" / "static" / "js"
SWEEP_SECONDS = 30


def create_app(registry=None, sweep_seconds=SWEEP_SECONDS):
    from .routers import (analysis, comparison, config, datasets, de, deconvolution, embedding, events, layers, load,
                          metadata, plugins, regions, summary, tiles, uploads)

    registry = SessionRegistry() if registry is None else registry

    @contextlib.asynccontextmanager
    async def lifespan(app):
        async def sweep():
            while True:
                await asyncio.sleep(sweep_seconds)
                registry.sweep()
        task = asyncio.create_task(sweep())
        try:
            yield
        finally:
            task.cancel()

    app = FastAPI(title="VisiumVisualizer API", lifespan=lifespan,
                  root_path=os.environ.get("VHD_ROOT_PATH", ""))
    app.state.sessions = registry

    @app.exception_handler(TooManySessions)
    async def busy(request, exc):
        return JSONResponse({"detail": "The server is busy; please try again in a few minutes."},
                            status_code=503, headers={"Retry-After": "120"})

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
            return FileResponse(target)
        return FileResponse(DIST / "index.html")

    return app


app = create_app()
