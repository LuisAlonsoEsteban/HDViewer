"""Every update the server sends must only reference models the browser knows.

BokehJS drops a patch that references an unknown model ("reference ... isn't
known"); new cluster widgets were lost that way."""
import json
import queue
import time
import unittest

from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import Button, CheckboxGroup, ImportedStyleSheet
from bokeh.protocol import Protocol

from hd_viewer.presentation import configure_document
from tests.test_clustering_markers import make_layer


class ServerLikeDoc(Document):
    """Next-tick callbacks run later, after the current change was sent."""
    def __init__(self):
        super().__init__()
        self.ticks = queue.Queue()

    def add_next_tick_callback(self, callback):
        self.ticks.put(callback)

    def add_periodic_callback(self, callback, period_ms):
        return callback

    def remove_periodic_callback(self, callback):
        pass

    def pump(self, seconds):
        end = time.monotonic() + seconds
        while time.monotonic() < end:
            try:
                self.ticks.get(timeout=.05)()
            except queue.Empty:
                pass


def _walk(value, defined, refs):
    if isinstance(value, dict):
        if value.get("type") == "bytes":          # binary buffer sent alongside
            return
        if "id" in value and len(value) == 1:
            refs.add(value["id"])
        elif "id" in value and value.get("type") == "object":
            defined.add(value["id"])
        for item in value.values():
            _walk(item, defined, refs)
    elif isinstance(value, (list, tuple)):       # to_json() may hold tuples
        for item in value:
            _walk(item, defined, refs)


class Browser:
    """Tracks the models a client would know from the patches it receives."""
    def __init__(self, doc):
        self.known, self.unknown = set(), []
        _walk(doc.to_json(), self.known, set())

    def _document_patched(self, event):
        content = json.loads(Protocol().create("PATCH-DOC", [event]).content_json)
        defined, refs = set(), set()
        _walk(content, defined, refs)
        missing = refs - self.known - defined
        if missing:
            self.unknown.append((type(event).__name__, getattr(event, "attr", None), sorted(missing)))
        self.known |= defined


class PatchOrderingTests(unittest.TestCase):
    def test_new_cluster_widgets_reach_the_browser_and_are_styled(self):
        layer, doc = make_layer(), ServerLikeDoc()
        configure_document(doc)                     # app code registers first, as in main.py
        root = layer.get_widgets_layout(doc)
        doc.add_root(root)
        doc.add_root(layer.figure)
        doc.pump(.5)
        browser = Browser(doc)
        doc.callbacks.on_change_dispatch_to(browser)   # then the server session

        button = next(b for b in root.select({"type": Button}) if b.label == "Run clustering")
        button._trigger_event(ButtonClick(button))
        cm = layer.clustering_manager
        deadline = time.monotonic() + 90
        while time.monotonic() < deadline and not any(
                isinstance(w, CheckboxGroup) for w in cm._cluster_layout.children):
            doc.pump(.2)
        doc.pump(1)

        boxes = [w for w in cm._cluster_layout.children if isinstance(w, CheckboxGroup)]
        self.assertTrue(boxes)
        self.assertEqual(browser.unknown, [])
        for box in boxes:                           # still styled, one tick later
            self.assertTrue(any(isinstance(s, ImportedStyleSheet) for s in box.stylesheets))


if __name__ == "__main__":
    unittest.main()
