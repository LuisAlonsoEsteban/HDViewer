"""Headless Chromium (Playwright) helpers shared by the browser-based tools.

A view counts as *rendered* when the page has stopped changing:
  * no network request (viewport data, image tiles, Zarr chunks) is pending, and
  * Chromium has composited no visually different frame,
for `quiet` seconds. Frames come from the DevTools screencast, which emits a
frame (with its compositor timestamp) whenever the page repaints; a spinner or
progressive tile loading keeps producing frames, so they count as not done.
The reported time runs from the action to the later of the last request's end
and the last visual change (plus two animation frames when nothing changed), so
the quiet window itself is not counted.
"""
import hashlib
import re
import time

CHROMIUM_ARGS = [
    "--enable-precise-memory-info",       # exact performance.memory values
    "--ignore-gpu-blocklist",
    "--use-angle=swiftshader",            # WebGL without a GPU (deck.gl, Bokeh WebGL)
    "--enable-unsafe-swiftshader",
]

RAF2 = "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(() => r(performance.now()))))"

FRAME_RECORDER = """() => {
  if (window.__bench) return;
  window.__bench = {frames: [], last: null, on: false};
  const tick = (t) => {
    const b = window.__bench;
    if (b.on && b.last !== null) b.frames.push(t - b.last);
    b.last = t;
    requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
}"""


class BrowserSession:
    def __init__(self, viewport=(1400, 900), ignore=(r"/api/events", r"/sockjs", r"hot-update"), headless=True):
        self.viewport, self.headless = viewport, headless
        self.ignore = [re.compile(p) for p in ignore]
        self.pending = set()
        self.last_activity = time.perf_counter()
        self.requests = 0
        self.bytes = 0
        self.failed = []

    def __enter__(self):
        from playwright.sync_api import sync_playwright
        self._pw = sync_playwright().start()
        self.browser = self._pw.chromium.launch(headless=self.headless, args=CHROMIUM_ARGS)
        self.context = self.browser.new_context(viewport=dict(width=self.viewport[0], height=self.viewport[1]))
        self.page = self.context.new_page()
        self.page.set_default_timeout(600_000)
        self.page.on("request", self._on_request)
        self.page.on("requestfinished", self._on_finished)
        self.page.on("requestfailed", self._on_failed)
        self.page.on("pageerror", lambda exc: self.failed.append(f"pageerror: {exc}"))
        self.last_visual = 0.0          # wall time (epoch s) of the last visually new frame
        self._last_digest = None
        self._cdp = self.context.new_cdp_session(self.page)
        self._cdp.on("Page.screencastFrame", self._on_frame)
        self._cdp.send("Page.startScreencast", dict(format="jpeg", quality=40, maxWidth=self.viewport[0],
                                                    maxHeight=self.viewport[1], everyNthFrame=1))
        return self

    def _on_frame(self, frame):
        try:
            self._cdp.send("Page.screencastFrameAck", dict(sessionId=frame["sessionId"]))
        except Exception:
            pass
        digest = hashlib.md5(frame["data"].encode()).digest()
        if digest != self._last_digest:
            self._last_digest = digest
            self.last_visual = max(self.last_visual, float(frame["metadata"]["timestamp"]))

    def __exit__(self, *exc):
        try:
            try:
                self._cdp.send("Page.stopScreencast")
            except Exception:
                pass
            self.browser.close()
        finally:
            self._pw.stop()

    # -- network tracking ---------------------------------------------------------
    def _tracked(self, request):
        return not any(p.search(request.url) for p in self.ignore)

    def _on_request(self, request):
        if self._tracked(request):
            self.pending.add(request)
            self.requests += 1
            self.last_activity = time.perf_counter()

    def _on_finished(self, request):
        if request in self.pending:
            self.pending.discard(request)
            self.last_activity = time.perf_counter()
            try:
                sizes = request.sizes()
                self.bytes += sizes.get("responseBodySize", 0) + sizes.get("responseHeadersSize", 0)
            except Exception:
                pass

    def _on_failed(self, request):
        if request in self.pending:
            self.pending.discard(request)
            self.last_activity = time.perf_counter()
            self.failed.append(f"{request.url}: {request.failure}")

    # -- timing ---------------------------------------------------------------------
    def settle(self, started, quiet=0.4, timeout=600.0, min_wait=0.15):
        """Seconds from `started` (perf_counter) until the view is rendered (see module doc)."""
        deadline = started + timeout
        offset = time.time() - time.perf_counter()      # screencast timestamps are wall-clock
        self.page.wait_for_timeout(min_wait * 1000)       # let the action issue its requests

        def last_change():
            return max(self.last_activity, self.last_visual - offset, started)

        while True:
            now = time.perf_counter()
            if not self.pending and now - last_change() >= quiet:
                break
            if now > deadline:
                raise TimeoutError(f"view not settled ({len(self.pending)} requests pending)")
            self.page.wait_for_timeout(20)
        end = last_change()
        if end > started:
            return end - started
        before = time.perf_counter()        # nothing changed: the cost of the next paint
        self.page.evaluate(RAF2)
        return time.perf_counter() - before

    def timed(self, action, **settle_kwargs):
        """Run `action()` and return (seconds until rendered, requests, bytes)."""
        requests0, bytes0 = self.requests, self.bytes
        started = time.perf_counter()
        action()
        seconds = self.settle(started, **settle_kwargs)
        return seconds, self.requests - requests0, self.bytes - bytes0

    # -- frames and memory ---------------------------------------------------------
    def frames_start(self):
        self.page.evaluate(FRAME_RECORDER)
        self.page.evaluate("() => { window.__bench.frames = []; window.__bench.last = null; window.__bench.on = true; }")

    def frames_stop(self):
        frames = self.page.evaluate("() => { window.__bench.on = false; return window.__bench.frames; }")
        if not frames:
            return {}
        frames = sorted(frames)
        n = len(frames)
        return dict(frames=n, frame_ms_median=frames[n // 2], frame_ms_p95=frames[min(n - 1, int(n * 0.95))],
                    frame_ms_max=frames[-1], long_frames_over_50ms=sum(f > 50 for f in frames))

    def js_heap(self):
        try:
            return self.page.evaluate(
                "() => performance.memory ? {used: performance.memory.usedJSHeapSize,"
                " total: performance.memory.totalJSHeapSize} : null")
        except Exception:
            return None


def wait_http(url, timeout=600.0):
    """Block until `url` answers (a server started by an adapter)."""
    import urllib.request
    deadline = time.monotonic() + timeout
    while True:
        try:
            with urllib.request.urlopen(url, timeout=5) as response:
                if response.status < 500:
                    return
        except Exception:
            if time.monotonic() > deadline:
                raise TimeoutError(f"{url} did not start")
            time.sleep(0.25)


def free_port():
    import socket
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]
