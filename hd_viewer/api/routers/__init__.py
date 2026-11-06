from fastapi import Request, Response


def current_session(request: Request, response: Response):
    """FastAPI dependency: the caller's session (created on first use)."""
    return request.app.state.sessions.resolve(request, response)
