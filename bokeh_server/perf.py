"""Timings of the map's loading pipeline, printed when VISIUM_PERF=1.

Each viewport load prints one line, e.g.

    [perf] Bins agg:4 +3120 rows (24946 shown, 15608 vertices): wait 2 ms, fetch 310 ms, prepare 12 ms, stream 85 ms

wait: queued before a worker thread took it (background loads only)
fetch: Sample.fetch_data (geometry, values; runs in the worker thread)
prepare: colours, alphas and filters of the new rows (Bokeh thread)
stream: handing the rows to Bokeh, which serialises them for the browser

What happens after that (transfer and drawing) is the browser's: measure it with the
Performance tab of its developer tools.
"""

import os

ENABLED = os.environ.get("VISIUM_PERF", "") not in ("", "0")

def log(label, n_new, n_shown, n_vertices, **timings_ms):
    if not ENABLED:
        return
    times = ", ".join(f"{name} {value:.0f} ms" for name, value in timings_ms.items() if value is not None)
    print(f"[perf] {label} +{n_new} rows ({n_shown} shown, {n_vertices} vertices): {times}", flush=True)
