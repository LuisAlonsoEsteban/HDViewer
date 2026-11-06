"""Server-to-browser events of a session (load progress, results, heartbeat).

The socket also ties the session's lifetime to the open page: when the last
socket closes, the session ends after a grace period (see api.sessions).
"""
import asyncio
import json
import time

from fastapi import APIRouter, WebSocket, WebSocketDisconnect

from ..sessions import COOKIE

router = APIRouter()
HEARTBEAT_SECONDS = 15


@router.websocket("/events")
async def events(websocket: WebSocket):
    # The session (and its cookie) comes from an earlier HTTP request, e.g. GET /api/load
    session = websocket.app.state.sessions.get(websocket.cookies.get(COOKIE, ""))
    if session is None:
        await websocket.close(code=4401)
        return
    await websocket.accept()
    queue = asyncio.Queue()
    session.loop = asyncio.get_running_loop()
    session.subscribers.add(queue)
    session.disconnected_at = None
    session.touch()
    await websocket.send_text(json.dumps(dict(type="status", load=session.load)))

    async def send():
        while True:
            try:
                event = await asyncio.wait_for(queue.get(), HEARTBEAT_SECONDS)
            except asyncio.TimeoutError:
                event = dict(type="ping")
            await websocket.send_text(json.dumps(event))

    async def receive():
        from .load import cancel_load
        while True:
            message = json.loads(await websocket.receive_text() or "{}")
            session.touch()
            if message.get("type") == "cancel":
                cancel_load(session)

    tasks = [asyncio.create_task(send()), asyncio.create_task(receive())]
    try:
        await asyncio.wait(tasks, return_when=asyncio.FIRST_EXCEPTION)
    except WebSocketDisconnect:
        pass
    finally:
        for task in tasks:
            task.cancel()
        session.subscribers.discard(queue)
        if not session.subscribers:
            session.disconnected_at = time.monotonic()
