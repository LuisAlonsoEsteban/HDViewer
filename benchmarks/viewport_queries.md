# Indexed viewport queries

Native bin, cell, and spot cropping now queries the existing GeoPandas spatial
index with the same exact intersection predicate previously used by `.cx`.
Selected positions are sorted to preserve original row order and duplicate labels.

Aggregated bins use a separate compact index per resolution: one point per
aggregate block and a positional map back to all native members. Queries select
nearby blocks, expand their complete membership, and apply the original bounds
inequalities to those candidates. A small conservative expansion protects
floating-point boundary contacts. Aggregation values and category tie ordering
remain identical. The block index is reused and rebuilt when its cached grid
changes; replacing the geometry frame clears the per-resolution indexes.

## Repeat-query results

Times are milliseconds, memory is MiB, and entries show before → after.
Lookup means selecting/copying candidate geometry rows (and the aggregate-grid
rows for aggregated bins). Complete fetch additionally prepares coordinates,
observation names, gene values, and aggregation. Query-result cache hits are excluded.

| Scenario | Native rows | Lookup time | Complete fetch time | Peak temporary lookup allocation | Peak temporary fetch allocation |
|---|---:|---:|---:|---:|---:|
| bins_native_256 | 65,536 | 2.70 → 0.61 | 3.96 → 1.92 | 0.15 → 0.16 | 0.21 → 0.21 |
| bins_native_768 | 589,824 | 20.85 → 0.62 | 22.64 → 2.00 | 0.65 → 0.16 | 0.65 → 0.21 |
| cells_native_768 | 589,824 | 23.35 → 0.59 | 25.67 → 2.84 | 0.65 → 0.15 | 0.99 → 0.99 |
| bins_aggregate4_768 | 589,824 | 3.66 → 1.11 | 9.76 → 7.11 | 5.63 → 0.32 | 5.63 → 0.56 |
| bins_rotated_aggregate4_768 | 589,824 | 3.45 → 1.10 | 9.35 → 7.06 | 5.63 → 0.34 | 5.63 → 0.62 |
| bins_aggregate4_wide_768 | 589,824 | 5.33 → 2.80 | 24.45 → 21.22 | 5.63 → 4.15 | 8.14 → 7.58 |

The native-bin view selects 1,764 rows at both sample sizes. Increasing the
sample from 65,536 to 589,824 rows increases the original lookup from 2.70 to
20.85 ms; indexed lookup stays approximately 0.61–0.62 ms.

The small aggregated views select complete blocks, including native bins
outside the viewport. Full fetch improvements are smaller than lookup gains
because value preparation and groupby aggregation still take time.

## Construction and persistent-memory tradeoff

Indexes make repeated queries faster at the cost of first-use construction and
retained memory. They are built lazily, once per native geometry array or
aggregation grid. The native index is shared with other GeoPandas spatial-index
users, including ROI membership queries.

| Scenario | First lookup time | Peak temporary first-lookup allocation | Observed RSS growth during first lookup |
|---|---:|---:|---:|
| bins_native_256 | 4.53 → 14.69 | 0.16 → 1.50 | 0.00 → 3.44 |
| bins_native_768 | 22.18 → 135.76 | 0.66 → 13.50 | 0.32 → 30.00 |
| cells_native_768 | 53.34 → 217.32 | 0.65 → 13.50 | 0.19 → 51.61 |
| bins_aggregate4_768 | 6.01 → 24.77 | 5.63 → 17.73 | 0.55 → 0.81 |
| bins_rotated_aggregate4_768 | 5.87 → 25.72 | 5.63 → 17.73 | 0.00 → 0.72 |
| bins_aggregate4_wide_768 | 7.67 → 26.70 | 5.63 → 17.73 | 0.00 → 1.06 |

RSS growth is measured within each process before and after the first query
and garbage collection. It is not the index object size: allocator reuse can
hide retained allocations in RSS. For the 589,824-row aggregate cases, the
retained positional map and block offsets alone occupy about 4.78 MiB per
resolution, plus the point geometries and spatial tree. Absolute process RSS
and process high-water marks are retained in the raw results; imports and
dataset construction make those unsuitable for attributing the index cost alone.

Small native views use slightly more temporary lookup memory after indexing;
cell full-fetch allocation peaks are essentially unchanged because coordinate
extraction dominates that peak. The benefit is lower repeated-query work and
fewer sample-sized temporary arrays, rather than lower total memory in every case.

## Reproduce

Run sequentially from the repository root in the same environment:

```bash
MPLCONFIGDIR=/tmp/visium-mpl .venv/bin/python benchmarks/viewport_queries.py --baseline-ref 3fcb959 --output benchmarks/viewport_queries_before.json
MPLCONFIGDIR=/tmp/visium-mpl .venv/bin/python benchmarks/viewport_queries.py --output benchmarks/viewport_queries_after.json
```

The baseline backend is loaded read-only from Git. No working files are replaced.
The benchmark uses actual Sample/Bins/Cells code, synthetic geometry and AnnData,
one numeric gene, and eight moving viewports per process. Standard windows span
80 plot units per axis, wide windows 320; native bins are 2 units wide. Rotated
aggregate grids use a 27-degree angle. Full sample gene, centroid, observation
mapping, and aggregation-grid caches are prepared equally before timing.
The additional viewport index is not prepared in advance.

Each scenario has three fresh timing processes. Reported warm timings are the
median of each process’s eight-view median. Memory is measured in a separate
fresh process with tracemalloc; warm peaks are the maximum across its eight
views. Tracemalloc covers tracked Python/NumPy allocations and does not capture
all GEOS/native allocations. RSS includes those native allocations.

All six scenarios, all timing repetitions, and allocation runs produce identical
ordered fetch-result hashes and selected-row counts before and after. These are
backend measurements; Bokeh source updates, network transfer, and browser rendering
are excluded. Dataset creation and common-cache preparation are also excluded.

Validation: 44 tests pass: 13 indexed-query tests, 15 source-update tests, six
parameter-form tests, and ten existing cell/viewport/detail-level tests. The
latter’s legacy imports were mapped to the current packages in memory.

Environment: Python 3.11.13; numpy 2.4.6, pandas 2.3.3, geopandas 1.1.4, shapely 2.1.2, anndata 0.12.19.
