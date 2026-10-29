"""Lightweight progress reporting for long-running computations.

Computation code calls `report("step...")`; the caller decides where messages go by wrapping the
work in `reporting(callback)`. The callback is thread-local, so concurrent jobs running in their own
threads (e.g. from the Bokeh server) each log to their own sink. Without a callback, messages are
ignored.
"""

import threading
from contextlib import contextmanager

_local = threading.local()


def report(message):
    callback = getattr(_local, "callback", None)
    if callback is not None:
        callback(message)


@contextmanager
def reporting(callback):
    previous = getattr(_local, "callback", None)
    _local.callback = callback
    try:
        yield
    finally:
        _local.callback = previous
