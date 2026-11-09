"""Server half of chunked file uploads, independent of the transport.

A file arrives as numbered raw-byte chunks after a "start" that announces its
name and size. The name is reduced to its final component and checked against
the allowed extensions, the size against the limit and the free disk space.
Chunks are appended to a temporary file (only one chunk is in memory) and the
complete file is renamed to ``<directory>/<key>/<name>``.
"""
import os
import shutil
from pathlib import Path

from .. import settings


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


class UploadField:
    """One upload slot (e.g. the expression matrix). `directory` is a callable
    returning the session folder. Methods return an error message or None."""

    EMPTY, UPLOADING, COMPLETE, ERROR = "empty", "uploading", "complete", "error"

    def __init__(self, key, suffixes, directory, max_bytes=None, chunk_bytes=None):
        self.key = key
        self.suffixes = tuple(suffix.lower() for suffix in suffixes)
        self.accept = ",".join(self.suffixes)
        self._directory = directory
        self.max_bytes = settings.UPLOAD_MAX_BYTES if max_bytes is None else max_bytes
        self.chunk_bytes = settings.UPLOAD_CHUNK_BYTES if chunk_bytes is None else chunk_bytes
        self._reset()

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
    def upload_id(self):
        return self._id

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
        return message

    def clear(self):
        """Cancel an upload in progress or delete the received file."""
        self._discard_files()
        self._reset()

    def begin(self, upload_id, name, size):
        """Start a new file (replacing any previous one)."""
        self._discard_files()
        self._reset()
        self._id = upload_id
        try:
            self.filename = safe_upload_name(name)
        except ValueError as exc:
            self.filename = str(name or "")
            return self._fail(str(exc))
        if not self.filename.lower().endswith(self.suffixes):
            return self._fail(f"Invalid file: {self.filename}. Expected {' or '.join(self.suffixes)}.")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            return self._fail(f"Invalid file: {self.filename} is empty.")
        if size > self.max_bytes:
            return self._fail(f"File too large: {self.filename} is {format_bytes(size)}; "
                              f"the limit is {format_bytes(self.max_bytes)}.")
        directory = Path(self._directory())
        try:
            (directory / self.key).mkdir(parents=True, exist_ok=True)
            if shutil.disk_usage(directory).free < size:
                return self._fail(f"Not enough disk space on the server for {self.filename}.")
            self._temp = directory / f".{self.key}.part"
            self._temp.write_bytes(b"")
        except OSError as exc:
            return self._fail(f"Cannot store {self.filename}: {exc.strerror or exc}")
        self.size = size
        self.status = self.UPLOADING
        return None

    def is_current(self, upload_id):
        return self.status == self.UPLOADING and upload_id == self._id

    def write(self, index, chunk):
        """Append chunk number `index` (bytes) of the current upload."""
        if index != self._next_index:
            return self._fail(f"Upload of {self.filename} failed: chunk {index} arrived "
                              f"out of order (expected {self._next_index}). Choose the file again.")
        if not chunk or len(chunk) > self.chunk_bytes:
            return self._fail(f"Upload of {self.filename} failed: invalid chunk size.")
        if self.received + len(chunk) > self.size:
            return self._fail(f"Upload of {self.filename} failed: more data than the "
                              f"announced {format_bytes(self.size)}.")
        try:
            with open(self._temp, "ab") as file:
                if file.tell() != self.received:
                    raise OSError("partial file changed on disk")
                file.write(chunk)
        except OSError as exc:
            return self._fail(f"Cannot store {self.filename}: {exc.strerror or exc}")
        self.received += len(chunk)
        self._next_index += 1
        if self.received == self.size:
            final = self._temp.parent / self.key / self.filename
            try:
                os.replace(self._temp, final)
            except OSError as exc:
                return self._fail(f"Cannot store {self.filename}: {exc.strerror or exc}")
            self._temp, self.path, self.status = None, final, self.COMPLETE
        return None

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
