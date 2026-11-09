"""Per-session server folders and their cleanup (shared by every front end)."""
import hashlib
import shutil
import sys

from .. import settings
from ..plugins.jobs import JOBS


def session_upload_dir(session_id):
    """Folder of a session's uploads and results. Session ids may come from the
    client, so never use one as a path: a hash cannot escape UPLOAD_DIR."""
    return settings.UPLOAD_DIR / hashlib.sha256(session_id.encode()).hexdigest()[:32]


def cleanup_session(session_id):
    """Stop the session's segmentation/deconvolution runs (unless another session
    awaits the same one), then remove the folder they write to."""
    cancelled = JOBS.cancel_owner(session_id)
    if cancelled:
        print(f"[session] Cancelled {cancelled} analysis job(s) of a closed session")
    upload_dir = session_upload_dir(session_id)
    if upload_dir.exists():
        try:
            shutil.rmtree(upload_dir)
            print(f"[session] Removed upload directory {upload_dir}")
        except Exception as exc:
            print(f"[session] Failed to remove upload directory {upload_dir}: {exc}", file=sys.stderr)
