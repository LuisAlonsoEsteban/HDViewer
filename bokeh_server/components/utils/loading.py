import html
import threading
import time
import traceback

from bokeh.models import Div

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from src.utils.progress import reporting

_SPINNER_CSS = """
<style>
@keyframes loading-logger-spin { to { transform: rotate(360deg); } }
.loading-logger-spinner {
    display: inline-block; width: 10px; height: 10px; margin-right: 6px;
    border: 2px solid #ccc; border-top-color: #3182bd; border-radius: 50%;
    animation: loading-logger-spin 0.8s linear infinite; vertical-align: middle;
}
.loading-logger-log { color: #888; font-size: 11px; margin: 2px 0 0 18px; }
</style>
"""


class LoadingLogger:
    """Status Div showing a spinner, the current step, elapsed time and a short log while a
    background job runs.

    `run(title, work, on_done)` executes `work()` in a thread; messages sent through
    `src.utils.progress.report` from inside it are shown live. `on_done(result)` runs on the
    document thread afterwards, so it may update Bokeh models.
    """

    MAX_LOG_LINES = 6

    def __init__(self, doc, text="", font_size=12, **div_kwargs):
        self.doc = doc
        self.font_size = font_size
        self.div = Div(text=text, **div_kwargs)
        self.running = False
        self._title = ""
        self._steps = []
        self._start = None
        self._periodic = None
        self._lock = threading.Lock()

    # ------------------------------------------------------------------
    # Static messages
    # ------------------------------------------------------------------

    def set_message(self, message, color="gray"):
        self.div.text = f"<p style='color:{color};font-size:{self.font_size}px;margin:0'>{message}</p>"

    def success(self, message):
        self.set_message(f"{message} <span style='color:#888'>({self._elapsed():.1f}s)</span>", color="green")

    def error(self, message):
        self.set_message(f"Error: {html.escape(str(message))}", color="red")

    # ------------------------------------------------------------------
    # Running a job
    # ------------------------------------------------------------------

    def run(self, title, work, on_done=None, on_error=None):
        """Run `work()` in a background thread while showing progress.

        on_done(result) / on_error(exc) are called on the document thread. Without `on_error`, the
        exception message is shown in the Div.
        """
        if self.running:
            return
        self.running = True
        self._title = title
        self._steps = []
        self._start = time.monotonic()
        self._render()
        self._periodic = self.doc.add_periodic_callback(self._render, 500)

        def _target():
            try:
                with reporting(self.log):
                    result = work()
            except Exception as exc:
                print(traceback.format_exc())
                self.doc.add_next_tick_callback(lambda exc=exc: self._finish(on_error, exc, failed=True))
            else:
                self.doc.add_next_tick_callback(lambda: self._finish(on_done, result, failed=False))

        threading.Thread(target=_target, daemon=True).start()

    def log(self, message):
        """Thread-safe: record a step (also printed to the server console)."""
        print(f"[{self._title}] {message}")
        with self._lock:
            self._steps.append((self._elapsed(), str(message)))

    def _finish(self, callback, value, failed):
        self.running = False
        if self._periodic is not None:
            try:
                self.doc.remove_periodic_callback(self._periodic)
            except ValueError:
                pass
            self._periodic = None

        if failed:
            if callback is not None:
                callback(value)
            else:
                self.error(value)
        elif callback is not None:
            callback(value)
        else:
            self.success("Done")

    def _elapsed(self):
        return 0.0 if self._start is None else time.monotonic() - self._start

    def _render(self):
        if not self.running:
            return
        with self._lock:
            steps = list(self._steps)

        current = html.escape(steps[-1][1]) if steps else "Starting..."
        previous = "".join(
            f"<div class='loading-logger-log'>{t:5.1f}s &nbsp;{html.escape(msg)}</div>"
            for t, msg in steps[-self.MAX_LOG_LINES:-1]
        )
        self.div.text = (
            f"{_SPINNER_CSS}"
            f"<div style='font-size:{self.font_size}px'>"
            f"{previous}"
            f"<div><span class='loading-logger-spinner'></span>"
            f"<b>{html.escape(self._title)}</b> &middot; {current} "
            f"<span style='color:#888'>({self._elapsed():.0f}s)</span></div>"
            f"</div>"
        )
