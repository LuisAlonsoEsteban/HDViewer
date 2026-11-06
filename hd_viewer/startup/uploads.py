"""Chunked file uploads over the Bokeh websocket.

Bokeh's FileInput sends a whole file, base64-encoded, in a single websocket
message, which the default 20 MB message limit rejects for real Visium HD data
(and which would otherwise stay in memory for the whole session). Instead the
browser (static/js/chunked_upload.js) reads the chosen file in slices and sends
one slice at a time through a hidden model property; the server appends each
slice to a temporary file and acknowledges it before the next one is sent.

Protocol (JSON strings on ``UploadChannel``):
  browser -> ``request``:  {"action": "start", "id", "name", "size"}
                           {"action": "chunk", "id", "index", "data"}  (base64)
  server  -> ``response``: {"seq", "id", "index", "ok", "done", "error"}
``index`` is -1 for the acknowledgement of "start".
"""

import base64
import binascii
import json
from pathlib import Path

from bokeh.core.properties import String
from bokeh.layouts import row
from bokeh.model import DataModel
from bokeh.models import Button, CustomJS

from ..core.uploads import UploadField, format_bytes, safe_upload_name  # noqa: F401  (re-exported)

_JS_PATH = Path(__file__).resolve().parent.parent / "static" / "js" / "chunked_upload.js"
UPLOAD_JS = _JS_PATH.read_text(encoding="utf-8")


class UploadChannel(DataModel):
    """Hidden model carrying one upload message in each direction."""
    request = String(default="")
    response = String(default="")


class ChunkedUpload(UploadField):
    """One upload field: picker buttons plus the server half of the protocol.

    The assembled file lives at ``directory()/<key>/<sanitized name>``; only one
    chunk is ever held in memory.
    """

    def __init__(self, key, suffixes, directory, on_change=None, can_upload=lambda: True,
                 max_bytes=None, chunk_bytes=None):
        super().__init__(key, suffixes, directory, max_bytes=max_bytes, chunk_bytes=chunk_bytes)
        self._on_change = on_change
        self._can_upload = can_upload
        self._seq = 0

        self.channel = UploadChannel()
        self.choose_button = Button(label="Choose file…", button_type="default")
        self.remove_button = Button(label="Remove", button_type="default", visible=False)
        self.widget = row(self.choose_button, self.remove_button, sizing_mode="stretch_width")

        self.channel.on_change("request", self._handle_request)
        self.channel.js_on_change("response", CustomJS(
            code="cb_obj.__hd_upload?.on_response(cb_obj.response);"))
        self.remove_button.on_click(self.remove)

    def bind_js(self, load_button, data_source):
        self.choose_button.js_on_click(CustomJS(args=dict(
            channel=self.channel, accept=self.accept, chunk_bytes=self.chunk_bytes,
            load_button=load_button, source=data_source,
        ), code=UPLOAD_JS))

    # ----------------------------------------------------------------- state
    @property
    def disabled(self):
        return self.choose_button.disabled

    @disabled.setter
    def disabled(self, value):
        self.choose_button.disabled = value
        self.remove_button.disabled = value

    def _changed(self):
        self.remove_button.visible = self.status != self.EMPTY
        self.choose_button.label = "Choose file…" if self.status == self.EMPTY else "Choose another file…"
        self.remove_button.label = "Cancel upload" if self.status == self.UPLOADING else "Remove"
        if self._on_change is not None:
            self._on_change()

    def remove(self):
        """Cancel an upload in progress or delete the received file."""
        if not self._can_upload():
            return
        upload_id = self._id
        self.clear()
        if upload_id is not None:
            self._respond(upload_id, -1, ok=False, error="cancelled")
        self._changed()

    # -------------------------------------------------------------- protocol
    def _respond(self, upload_id, index, ok=True, done=False, error=None):
        self._seq += 1
        self.channel.response = json.dumps(dict(
            seq=self._seq, id=upload_id, index=index, ok=ok, done=done,
            error=error, received=self.received))

    def _handle_request(self, attr, old, new):
        if not new:
            return
        # Drop the chunk from the model as soon as possible.
        self.channel.request = ""
        try:
            message = json.loads(new)
            action, upload_id = message["action"], str(message["id"])
        except (ValueError, KeyError, TypeError):
            return
        if action == "start":
            self.start(upload_id, message.get("name"), message.get("size"))
        elif action == "chunk":
            self.receive_chunk(upload_id, message.get("index"), message.get("data"))

    def start(self, upload_id, name, size):
        if not self._can_upload():
            self._respond(upload_id, -1, ok=False, error="Uploads are locked while data is loading.")
            return
        # Starting a new file always replaces (and blocks loading on) the previous one.
        error = self.begin(upload_id, name, size)
        self._respond(upload_id, -1, ok=error is None, error=error)
        self._changed()

    def receive_chunk(self, upload_id, index, data):
        if not self.is_current(upload_id):
            # A late chunk from a cancelled or replaced upload: ignore it.
            self._respond(upload_id, index, ok=False, error="This upload is no longer active.")
            return
        if index == self._next_index:
            # 4 base64 characters encode 3 bytes; check before decoding.
            if not isinstance(data, str) or not data or len(data) > 4 * (self.chunk_bytes // 3 + 1):
                return self._abort(index, self._fail(f"Upload of {self.filename} failed: invalid chunk size."))
            try:
                chunk = base64.b64decode(data, validate=True)
            except (binascii.Error, ValueError):
                return self._abort(index, self._fail(f"Upload of {self.filename} failed: corrupted chunk."))
        else:
            chunk = b""                    # reported as out of order by write()
        error = self.write(index, chunk)
        del chunk
        if error is not None:
            return self._abort(index, error)
        self._respond(upload_id, index, done=self.complete)
        self._changed()

    def _abort(self, index, message):
        self._respond(self._failed_id, index, ok=False, error=message)
        self._changed()

    def _fail(self, message):
        self._failed_id = self._id
        return super()._fail(message)
