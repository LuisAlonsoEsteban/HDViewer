"""Binary frames carrying layer data to the browser.

Layout: a little-endian uint32 with the header length, the JSON header, then
numeric buffers, each aligned to 8 bytes so the browser can view them as typed
arrays without copying. A column is either {"dtype", "offset", "length"}
(binary) or {"values": [...]} (JSON: strings, polygon rings).

The header lists `ops` (as recorded by core.layer_state.RecordingSink):
{"op": "append" | "replace", "columns": {name: column}} or
{"op": "set", "name": ..., "column": column}.
"""
import json
import struct

import numpy as np

DTYPES = {np.dtype(np.float32): "float32", np.dtype(np.float64): "float64",
          np.dtype(np.uint8): "uint8", np.dtype(np.uint16): "uint16", np.dtype(np.uint32): "uint32",
          np.dtype(np.int8): "int8", np.dtype(np.int16): "int16", np.dtype(np.int32): "int32",
          np.dtype(np.bool_): "bool"}


def _json_values(values):
    out = []
    for value in values:
        if isinstance(value, np.ndarray):
            out.append(value.tolist())
        elif isinstance(value, (np.floating, float)):
            out.append(None if np.isnan(value) else float(value))
        elif isinstance(value, np.integer):
            out.append(int(value))
        else:
            out.append(value)
    return out


class FrameWriter:
    def __init__(self):
        self.buffers = []
        self.offset = 0

    def column(self, values):
        if isinstance(values, np.ndarray) and values.dtype in DTYPES and values.ndim == 1:
            data = np.ascontiguousarray(values).tobytes()
            entry = dict(dtype=DTYPES[values.dtype], offset=self.offset, length=int(values.size))
            pad = (-len(data)) % 8
            self.buffers.append(data + b"\0" * pad)
            self.offset += len(data) + pad
            return entry
        if isinstance(values, np.ndarray) and values.dtype.kind in "iu":
            return self.column(values.astype(np.int32))
        if isinstance(values, np.ndarray) and values.dtype.kind == "f":
            return self.column(values.astype(np.float32))
        return dict(values=_json_values(values))

    def ops(self, ops):
        encoded = []
        for op in ops:
            if op[0] in ("append", "replace"):
                encoded.append(dict(op=op[0], columns={name: self.column(values) for name, values in op[1].items()}))
            elif op[0] == "set":
                encoded.append(dict(op="set", name=op[1], column=self.column(op[2])))
        return encoded


def encode(header, ops=(), extra=None):
    """Frame bytes for `header` (a dict) plus recorded `ops` and named `extra` arrays."""
    writer = FrameWriter()
    header = dict(header, ops=writer.ops(ops))
    if extra:
        header["extra"] = {name: writer.column(values) for name, values in extra.items()}
    text = json.dumps(header, separators=(",", ":")).encode()
    pad = (-(4 + len(text))) % 8
    return struct.pack("<I", len(text) + pad) + text + b" " * pad + b"".join(writer.buffers)


def decode(frame):
    """(header, columns resolved to numpy arrays / lists): for tests and debugging."""
    (length,) = struct.unpack_from("<I", frame)
    header = json.loads(frame[4:4 + length])
    body = memoryview(frame)[4 + length:]
    lookup = {name: np.dtype(np.bool_ if name == "bool" else name) for name in DTYPES.values()}

    def resolve(column):
        if "values" in column:
            return column["values"]
        dtype = lookup[column["dtype"]]
        return np.frombuffer(body, dtype=dtype, count=column["length"], offset=column["offset"])

    ops = []
    for op in header["ops"]:
        if op["op"] == "set":
            ops.append(("set", op["name"], resolve(op["column"])))
        else:
            ops.append((op["op"], {name: resolve(column) for name, column in op["columns"].items()}))
    extra = {name: resolve(column) for name, column in header.get("extra", {}).items()}
    return header, ops, extra
