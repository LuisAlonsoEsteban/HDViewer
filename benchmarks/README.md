# Source update performance

For the subsequent spatial-index work, see [viewport query benchmarks](viewport_queries.md).

The changes retain numeric source columns as NumPy arrays, preserve loaded
geometry during compatible parameter changes, and keep non-hovered parameter
values on the server for slider updates. Pruning preserves numeric arrays and
polygon rings while keeping raw values, sequencing masks, and selections aligned.
Unused parameter values are dropped when parameters are removed; all raw values
are cleared when the source context changes.

Gene changes fetch values without coordinates and align them by row identifier.
A resolution change, in-flight load, or incomplete coverage of retained rows
falls back to the existing full reload path. The value-only sample fetch now
retains its observation index so these values can be aligned correctly.

## Reproduce

Run from the repository root, using the same Python environment for both:

```bash
.venv/bin/python benchmarks/source_updates.py --baseline-ref cb7e366 --output benchmarks/source_updates_before.json
.venv/bin/python benchmarks/source_updates.py --output benchmarks/source_updates_after.json
```

The baseline manager is loaded read-only from Git, without replacing working
files. Both runs use the same benchmark and parameter implementation.

## Method

- Synthetic input has 60,000 observations. The initial viewport loads 40,000;
  streaming adds the remaining 20,000. Pruning removes every fourth loaded row,
  leaving 30,000. Style, scale, and gene changes operate on 40,000 rows.
- Bin coordinates are numeric; cells have 13 coordinates per polygon ring.
  Every eleventh gene value is missing, exercising unsequenced styling.
- The source manager and Bokeh's actual `PATCH-DOC` serializer are used. The
  synthetic sample supplies already prepared values, isolating source-update
  overhead from disk reads, spatial queries, and aggregation.
- Each timing repetition uses a fresh process. Reported times are medians of
  three runs. The source operation and protocol serialization are timed
  separately; their combined median is computed from each run's sum.
- Source-change payloads are snapshotted when their events occur, including
  copies of mutable columns. This prevents a clear-source event from accidentally
  serializing rows appended later. Snapshot capture is included in operation
  time and allocation measurements.
- Peak traced allocations come from a separate run with `tracemalloc`, covering
  the source operation and event capture, excluding later serialization. This
  captures tracked Python/NumPy allocations, not every native allocation.
- Process RSS is also recorded. Its high-water mark includes imports, synthetic
  sample construction, and setup through the source operation. Reported process
  peaks are medians of the three timing runs; they exclude later serialization.
- Patch size includes uncompressed content JSON and binary buffer bodies. It
  excludes protocol headers, network framing, and compression.
- Canonical hashes of all final source values are compared across versions.
  JSON results retain individual runs, row counts, fetch counts, serialization
  times, memory measurements, and package versions.

These are source-path benchmarks, not end-to-end app measurements. They exclude
image loading, analysis, real sample fetch costs, renderer reconstruction,
WebSocket transfer, browser memory, and rendering time. No runtime timing logger
is enabled by these changes.

## Results

The table shows before → after. Times are milliseconds; memory and patch sizes
are MiB. “Combined” includes the source operation and protocol serialization.

| Scenario | Source time | Combined time | Peak traced allocation | Process peak RSS | Patch size |
|---|---:|---:|---:|---:|---:|
| bins/initial | 130.2 → 92.8 | 308.0 → 134.2 | 17.5 → 11.5 | 161.9 → 154.8 | 1.85 → 2.36 |
| bins/stream | 70.3 → 54.5 | 158.3 → 76.4 | 11.5 → 9.8 | 165.8 → 157.9 | 0.93 → 1.18 |
| bins/prune | 58.3 → 63.2 | 192.9 → 94.4 | 4.8 → 4.8 | 161.4 → 154.3 | 1.38 → 1.77 |
| bins/style | 122.3 → 15.4 | 302.5 → 36.9 | 17.5 → 3.8 | 172.7 → 156.5 | 1.85 → 0.38 |
| bins/scale | 124.5 → 18.3 | 312.3 → 18.5 | 17.5 → 3.8 | 172.8 → 157.8 | 2.11 → 0.31 |
| bins/gene | 125.7 → 40.3 | 307.5 → 60.7 | 17.5 → 10.0 | 172.9 → 159.7 | 1.83 → 0.69 |
| cells/initial | 145.8 → 98.4 | 1227.0 → 1098.2 | 15.7 → 12.1 | 229.3 → 224.8 | 11.80 → 12.21 |
| cells/stream | 78.1 → 56.9 | 613.2 → 556.7 | 10.6 → 9.8 | 232.0 → 227.4 | 5.90 → 6.11 |
| cells/prune | 115.3 → 61.6 | 995.2 → 833.3 | 17.1 → 4.8 | 242.6 → 224.3 | 8.85 → 9.16 |
| cells/style | 136.0 → 15.4 | 1225.2 → 35.9 | 15.7 → 3.8 | 237.7 → 228.4 | 11.80 → 0.38 |
| cells/scale | 136.9 → 18.0 | 1223.9 → 18.2 | 15.7 → 3.8 | 237.4 → 226.7 | 12.07 → 0.31 |
| cells/gene | 138.3 → 42.9 | 1212.6 → 63.3 | 15.7 → 10.0 | 237.6 → 229.1 | 11.78 → 0.69 |

All 12 scenarios retain identical final source values and row counts.
Static style changes and non-hovered alpha scaling perform zero sample
fetches after the change; gene switches still fetch values once.

The largest gains come from avoiding full geometry reloads during style,
scale, and compatible gene changes. Polygon pruning also avoids the former
rectangular object-array expansion of equal-length rings.

There are measured tradeoffs: initial/streamed bin patches are about 28%
larger for these short numeric values because float64 binary buffers can
exceed their short JSON representations. Cell geometry patches grow about
3.5%. Precision is preserved. Serialization still becomes faster, and the
combined server times improve in every scenario. The isolated bin-pruning
operation becomes slightly slower, despite faster combined update/encoding.

Validation: 15 source-update tests, 6 parameter-form tests, and 10 existing
cell/viewport/detail-level tests. The latter use legacy module names; their
imports were mapped to current packages in memory for the regression run.
