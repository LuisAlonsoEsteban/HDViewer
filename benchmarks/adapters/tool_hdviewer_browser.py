"""HDViewer (ours): rendering in headless Chromium (run with the bench-sd env).

Starts the real server (`hdviewer_server.py` under the project's .venv), loads
the benchmark dataset through the loading screen like a user would, and times
each view until rendered (harness/browser.py). Views are set exactly by
writing the Bokeh map's ranges, which fires the same events as a pan/zoom.
Memory covers the whole tree: server + Chromium.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH.parent))

from benchmarks.harness.browser import BrowserSession, free_port, wait_http  # noqa: E402
from benchmarks.harness.worker import main  # noqa: E402

DATASET_TITLE = "Benchmark dataset"
CELLS_TITLE = "Cells · Cellpose"

SET_RANGES = """([x0, x1, y0, y1]) => {
  const plot = window.hdv.docs.get("map").get_model_by_name("map");
  plot.x_range.setv({start: x0, end: x1});
  plot.y_range.setv({start: y0, end: y1});
}"""


class Server:
    def __init__(self, spec, with_cells):
        self.spec = dict(spec, with_cells=with_cells, boxes_out=os.path.join(spec["workdir"], "map_boxes.json"))
        self.port = free_port()

    def __enter__(self):
        python = self.spec["hdviewer_python"]
        log = open(os.path.join(self.spec["workdir"], "hdviewer_server.log"), "w")
        self.proc = subprocess.Popen([python, str(BENCH / "adapters" / "hdviewer_server.py"),
                                      "--spec", json.dumps(self.spec), "--port", str(self.port)],
                                     stdout=log, stderr=subprocess.STDOUT, cwd=str(BENCH.parent))
        self.url = f"http://127.0.0.1:{self.port}"
        return self

    def wait(self):
        wait_http(self.url + "/api/datasets", timeout=1800)
        data = json.loads(Path(self.spec["boxes_out"]).read_text())
        return data["boxes"], data["seconds"]

    def __exit__(self, *exc):
        self.proc.terminate()
        try:
            self.proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            self.proc.kill()


def _load(b, url, with_cells):
    page = b.page
    # The cookie notice is fixed at the bottom of the page and would cover "Load data"
    page.add_init_script("try { sessionStorage.setItem('hdv:cookie-notice', 'seen') } catch (e) {}")
    page.goto(url + "/")
    page.get_by_role("button", name="Explore your data").click()
    select = page.get_by_role("combobox", name="Example dataset")
    for _ in range(20):                      # the dataset list arrives asynchronously
        select.click()
        page.get_by_role("option", name=DATASET_TITLE).click()
        if select.input_value() == DATASET_TITLE:
            break
        page.wait_for_timeout(250)
    page.get_by_label("Generate Visium v2 pseudospots (55 µm)").set_checked(False)
    if with_cells:
        page.get_by_label("Cellpose", exact=True).set_checked(True)
    with page.expect_response(lambda r: "/viewport" in r.url and r.ok, timeout=3_600_000):
        page.get_by_role("button", name="Load data").click()


def _colour(b, gene):
    page = b.page
    box = page.get_by_role("combobox", name="Gene / metadata")
    box.fill(gene)
    page.get_by_role("option", name=gene, exact=True).click()


def _session(spec, rec, modality):
    with_cells = modality == "cells"
    if with_cells:      # cells pool the 2 µm bins that fall within them
        spec = dict(spec, bin_size=spec["cells_bin_size"])
    with Server(spec, with_cells) as server:
        started = time.perf_counter()
        boxes, prep = server.wait()
        rec.add("server_start", time.perf_counter() - started - prep)
        with BrowserSession() as b:
            seconds, n, nbytes = b.timed(lambda: _load(b, server.url, with_cells), timeout=3600)
            rec.add("load_to_first_view", seconds, requests=n, bytes=nbytes)
            if modality == "e2e":
                rec.note(js_heap=b.js_heap())
                return
            page = b.page
            if modality in ("image", "cells"):
                page.get_by_role("button", name="Hide bins").click()
                b.settle(time.perf_counter())
            if modality == "cells":
                page.get_by_label(f"{CELLS_TITLE} layer").click()
            if modality != "image":
                seconds, n, nbytes = b.timed(lambda: _colour(b, spec["genes"][0]))
                rec.add("colour_by_gene", seconds, requests=n, bytes=nbytes)
            for name in [k for k in boxes if k != "pan"]:
                seconds, n, nbytes = b.timed(lambda: page.evaluate(SET_RANGES, boxes[name]))
                rec.add(f"view:{name}", seconds, requests=n, bytes=nbytes)
            page.evaluate(SET_RANGES, boxes["zoom_250um"])
            b.settle(time.perf_counter())
            b.frames_start()
            steps, total_n, total_bytes = [], 0, 0
            for box in boxes["pan"]:
                seconds, n, nbytes = b.timed(lambda: page.evaluate(SET_RANGES, box))
                steps.append(seconds)
                total_n += n
                total_bytes += nbytes
            frames = b.frames_stop()
            rec.add("pan", sum(steps), steps=len(steps), step_seconds=steps, requests=total_n, bytes=total_bytes, **frames)
            if modality != "image":
                seconds, n, nbytes = b.timed(lambda: _colour(b, spec["genes"][1]))
                rec.add("gene_switch", seconds, requests=n, bytes=nbytes)
            if spec.get("snapshots"):
                page.screenshot(path=os.path.join(spec["workdir"], f"hdviewer_{modality}.png"))
            rec.note(js_heap=b.js_heap(), page_errors=b.failed[:5])


TASKS = {
    "e2e.first_view": lambda spec, rec: _session(spec, rec, "e2e"),
    "image.render": lambda spec, rec: _session(spec, rec, "image"),
    "bins.render": lambda spec, rec: _session(spec, rec, "bins"),
    "cells.render": lambda spec, rec: _session(spec, rec, "cells"),
}

if __name__ == "__main__":
    main(TASKS, versions=("playwright",))
