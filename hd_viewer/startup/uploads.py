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
import os
import shutil
from pathlib import Path

from bokeh.core.properties import String
from bokeh.layouts import row
from bokeh.model import DataModel
from bokeh.models import Button, CustomJS

from .. import settings

_JS_PATH = Path(__file__).resolve().parent.parent / "static" / "js" / "chunked_upload.js"
UPLOAD_JS = _JS_PATH.read_text(encoding="utf-8")


class UploadChannel(DataModel):
    """Hidden model carrying one upload message in each direction."""
    request = String(default="")
    response = String(default="")


def safe_upload_name(name):
    """Client-supplied names are reduced to their final component."""
    name = Path(str(name).replace("\\", "/")).name
    if name in ("", ".", "..") or any(ord(char) < 32 for char in name):
        raise ValueError("Invalid upload filename")
    return name


def format_bytes(size):
    for unit in ("bytes", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "bytes" else f"{size:.1f} {unit}"
        size /= 1024


class ChunkedUpload:
    """One upload field: picker buttons plus the server half of the protocol.

    The assembled file lives at ``directory()/<key>/<sanitized name>``; only one
    chunk is ever held in memory.
    """

    EMPTY, UPLOADING, COMPLETE, ERROR = "empty", "uploading", "complete", "error"

    def __init__(self, key, suffixes, directory, on_change=None, can_upload=lambda: True,
                 max_bytes=None, chunk_bytes=None):
        self.key = key
        self.suffixes = tuple(suffix.lower() for suffix in suffixes)
        self.accept = ",".join(self.suffixes)
        self._directory = directory
        self._on_change = on_change
        self._can_upload = can_upload
        self.max_bytes = settings.UPLOAD_MAX_BYTES if max_bytes is None else max_bytes
        self.chunk_bytes = settings.UPLOAD_CHUNK_BYTES if chunk_bytes is None else chunk_bytes
        self._seq = 0
        self._reset()

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
    def _reset(self):
        self.status = self.EMPTY
        self.filename = ""
        self.size = 0
        self.received = 0
        self.path = None
        self.error = None
        self._id = None
        self._next_index = 0
        self._temp = None

    @property
    def complete(self):
        return self.status == self.COMPLETE

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

    def _discard_files(self):
        for path in (self._temp, self.path):
            if path is not None:
                try:
                    Path(path).unlink(missing_ok=True)
                except OSError:
                    pass

    def _fail(self, message):
        self._discard_files()
        name = self.filename
        self._reset()
        self.status, self.filename, self.error = self.ERROR, name, message

    def remove(self):
        """Cancel an upload in progress or delete the received file."""
        if not self._can_upload():
            return
        upload_id = self._id
        self._discard_files()
        self._reset()
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
        self._discard_files()
        self._reset()
        self._id = upload_id
        try:
            self.filename = safe_upload_name(name)
        except ValueError as exc:
            self.filename = str(name or "")
            return self._reject(upload_id, str(exc))
        if not self.filename.lower().endswith(self.suffixes):
            return self._reject(upload_id, f"Invalid file: {self.filename}. "
                                           f"Expected {' or '.join(self.suffixes)}.")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            return self._reject(upload_id, f"Invalid file: {self.filename} is empty.")
        if size > self.max_bytes:
            return self._reject(upload_id, f"File too large: {self.filename} is {format_bytes(size)}; "
                                           f"the limit is {format_bytes(self.max_bytes)}.")
        directory = Path(self._directory())
        try:
            (directory / self.key).mkdir(parents=True, exist_ok=True)
            if shutil.disk_usage(directory).free < size:
                return self._reject(upload_id, f"Not enough disk space on the server for {self.filename}.")
            self._temp = directory / f".{self.key}.part"
            self._temp.write_bytes(b"")
        except OSError as exc:
            return self._reject(upload_id, f"Cannot store {self.filename}: {exc.strerror or exc}")
        self.size = size
        self.status = self.UPLOADING
        self._respond(upload_id, -1)
        self._changed()

    def _reject(self, upload_id, message):
        self._fail(message)
        self._respond(upload_id, -1, ok=False, error=message)
        self._changed()

    def receive_chunk(self, upload_id, index, data):
        if self.status != self.UPLOADING or upload_id != self._id:
            # A late chunk from a cancelled or replaced upload: ignore it.
            self._respond(upload_id, index, ok=False, error="This upload is no longer active.")
            return
        if index != self._next_index:
            return self._abort(index, f"Upload of {self.filename} failed: chunk {index} arrived "
                                      f"out of order (expected {self._next_index}). Choose the file again.")
        # 4 base64 characters encode 3 bytes; check before decoding.
        if not isinstance(data, str) or not data or len(data) > 4 * (self.chunk_bytes // 3 + 1):
            return self._abort(index, f"Upload of {self.filename} failed: invalid chunk size.")
        try:
            chunk = base64.b64decode(data, validate=True)
        except (binascii.Error, ValueError):
            return self._abort(index, f"Upload of {self.filename} failed: corrupted chunk.")
        if self.received + len(chunk) > self.size:
            return self._abort(index, f"Upload of {self.filename} failed: more data than the "
                                      f"announced {format_bytes(self.size)}.")
        try:
            with open(self._temp, "ab") as file:
                if file.tell() != self.received:
                    raise OSError("partial file changed on disk")
                file.write(chunk)
        except OSError as exc:
            return self._abort(index, f"Cannot store {self.filename}: {exc.strerror or exc}")
        self.received += len(chunk)
        del chunk
        self._next_index += 1
        done = self.received == self.size
        if done:
            final = self._temp.parent / self.key / self.filename
            try:
                os.replace(self._temp, final)
            except OSError as exc:
                return self._abort(index, f"Cannot store {self.filename}: {exc.strerror or exc}")
            self._temp, self.path, self.status = None, final, self.COMPLETE
        self._respond(upload_id, index, done=done)
        self._changed()

    def _abort(self, index, message):
        upload_id = self._id
        self._fail(message)
        self._respond(upload_id, index, ok=False, error=message)
        self._changed()

    # ------------------------------------------------------------ messages
    def describe(self, required):
        if self.status == self.ERROR:
            return self.error
        if self.status == self.UPLOADING:
            percent = int(100 * self.received / self.size) if self.size else 0
            return (f"Uploading {self.filename}… {percent}% "
                    f"({format_bytes(self.received)} of {format_bytes(self.size)}).")
        if self.status == self.COMPLETE:
            return f"Uploaded: {self.filename} ({format_bytes(self.size)})"
        return "Required — waiting for a file." if required else "Optional — no file selected."
