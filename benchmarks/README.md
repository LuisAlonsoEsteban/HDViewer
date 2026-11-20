# Benchmarks: HDViewer vs napari, SpatialData, Vitessce and Loupe Browser

This suite measures **time** and **memory** on the same Visium HD data and the same views for five tools:

- HDViewer (this repository)
- [SpatialData](https://spatialdata.scverse.org) (spatialdata, spatialdata-io, spatialdata-plot)
- [napari](https://napari.org)
- [Vitessce](http://vitessce.io)
- [Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser)

It covers three kinds of data:

- the **image**: the full-resolution H&E;
- **bins**: 2, 8 and 16 µm;
- **cells**: segmented outlines with their pooled expression.

Tasks cover computation (reading, pyramids, aggregation, pooling, clustering, DE), data queries for a view, and rendering.

All four automated tools run on Linux, headless, on one node:

- Python tools run in-process.
- napari runs on an Xvfb X server with Mesa OpenGL.
- HDViewer's web interface and Vitessce run in headless Chromium, driven by Playwright.

Loupe Browser is a closed Windows/macOS desktop application, so it is measured by hand with a [protocol and a memory logger](adapters/loupe/PROTOCOL.md), and the report merges those rows.

## Tasks

Every tool implements what it can. A task without an equivalent is recorded as **n.a.** with the reason, so coverage is itself a result (`fig_outcomes_*`).

| Task | What is timed | HDViewer | SpatialData | napari | Vitessce |
|---|---|---|---|---|---|
| `image.read` | full-resolution image into memory | `VisiumHD.read_img` | dask-image → `.compute()` | `viewer.open` | n.a. |
| `image.pyramid` | multiscale image for viewing | align + WebP tile pyramid (`create_tiles`) | multiscale `Image2DModel` → Zarr | n.a. (uses SpatialData's) | n.a. (uses SpatialData's) |
| `image.query` | pixels/tiles one view needs | tile lookup + read | `bounding_box_query` + level read | n.a. | n.a. |
| `image.render` | views drawn: overview, 1 mm, 250 µm, 50 µm, 10-step pan | Chromium | spatialdata-plot (Agg) | `screenshot()` | Chromium |
| `bins.load` | Space Ranger outputs → tool's data | read + align + bin geometry + layer | `visium_hd` (+ Zarr write/read) | n.a. | conversion to a Vitessce store |
| `bins.query` | viewport → data for that view | `update_viewport` + binary frame | `bounding_box_query` + gene values | n.a. | n.a. |
| `bins.render` | gene colouring, views, pan, gene switch | Chromium | spatialdata-plot | Points layer | Chromium |
| `bins.aggregate` | 2 µm → 8/16 µm | multiscale aggregation | `aggregate`, `rasterize_bins` | n.a. | n.a. |
| `bins.cluster` | normalize, HVG, PCA, kNN, Leiden, UMAP | `ClusteringState` | scanpy | n.a. | n.a. |
| `bins.de` | Wilcoxon, 4 tissue quadrants vs rest | `DEAnalysis` | `rank_genes_groups` | n.a. | n.a. |
| `cells.load` | cell outlines (`cells.geojson`) | read + transform | `ShapesModel` | n.a. | label image + store |
| `cells.pool` | 2 µm bins summed per cell | `CellsPooler` | `spatialdata.aggregate` | n.a. | n.a. |
| `cells.query` / `cells.render` / `cells.cluster` / `cells.de` | as for bins | ✓ | ✓ | render only (Shapes) | render only (labels) |
| `e2e.first_view` | launch → first map with data | server + loading screen | read + first plot | read + first screenshot | conversion + page |

Each task records **phases**, for example `view:overview`, `pan`, `gene_switch`, `pca` or `leiden`. Phases named `setup:*` prepare inputs and are excluded from the results, for example building a Zarr cache the task starts from. Rendering tasks use the same boxes in full-resolution pixels for every tool (`harness/viewports.py`):

- the tissue overview;
- 1 mm, 250 µm and 50 µm windows on the tissue;
- ten pan steps of half a 250 µm window.

**When a view is rendered**

- **Browsers (HDViewer, Vitessce).** No request is pending and Chromium has composited no new frame for a short quiet period, observed through the DevTools screencast. The time runs to the last request or last visual change, not to the end of the quiet window.
- **napari.** `viewer.screenshot()` returns. It draws synchronously, including the multiscale tiles the view needs.
- **SpatialData.** `spatialdata-plot` has drawn the figure on an Agg canvas of the same size.

**Memory**

- Memory is the peak resident set of the **whole process tree**, sampled at 20 Hz. That covers the Python process, the HDViewer server, Chromium's processes and Xvfb clients.
- Per phase, the report also records the peak minus the value just before the phase.
- Chromium's shared pages are counted once per process, so browser figures are upper bounds.
- JS heap sizes are recorded as well.

**Fairness**

- Each run is a fresh process, so no state carries over between tasks or repeats.
- Each tool uses its documented workflow. One Space Ranger folder feeds every tool.
- Threads are pinned (`threads` in `config.yaml`).
- Analysis parameters are identical: min cells 3, target sum 10⁴, Seurat HVG, 30 PCs, k = 15, Leiden resolution 0.9.
- The first repeat is reported separately as "first", with OS caches still cold for that tool. The rest give the median and IQR.
- A per-task timeout (default 30 min) and a memory cap (default 100 GB) end a run with status `timeout` or `oom`.

**Known differences, stated in the paper**

- Pooling: HDViewer counts the bins inside a cell. `spatialdata.aggregate` counts bins that intersect it.
- HDViewer cells need 2 µm bins.
- Vitessce shows bins as a rasterized label image, as in Vitessce's own Visium HD example.
- napari shows bins as a Points layer, because a Shapes layer of millions of squares does not load.

## Data

- **Mouse Brain**: `data/Visium_HD_Mouse_Brain`. It has 6.3 M in-tissue 2 µm bins, 19,059 genes, a 23,947 × 18,872 px H&E, and 29,525 cellpose cells (`segmentation/cellpose/…/cells.geojson`).
- **Synthetic sweep**: `harness/synthetic.py` writes Space Ranger-format folders with 100k, 500k, 2M and 6M in-tissue 2 µm bins. Each also has 8/16 µm sums, an H&E-like image, cells and spatial expression domains. Folders are created on first use in `data/benchmark_synthetic/<size>`.

spatialdata-io's reader requires `<id>_feature_slice.h5` and `<id>_tissue_image.*`. `harness/datasets.py` builds a mirror of symlinks with those names plus a stub feature slice holding only the required `file_format` metadata, so nothing is written into the data folders.

## Setup

```bash
ENVS=/scratch/me/envs bash benchmarks/setup.sh
```

This creates `bench-sd` (SpatialData, Vitessce, Playwright), `bench-napari` (napari with PySide6) and `bench-xvfb` (X server for OpenGL). It also builds the Vitessce page (`adapters/vitessce_app`, with node from the `hdv-node` conda env). HDViewer runs in the project's `.venv`. Edit the paths in `config.yaml` if your `ENVS` differs.

## Run

```bash
.venv/bin/python benchmarks/run.py --list                                   # planned runs
.venv/bin/python benchmarks/run.py --datasets syn_100k --repeats 1          # quick check (~15 min)
.venv/bin/python benchmarks/run.py --datasets mouse_brain                   # real data, 5 repeats
.venv/bin/python benchmarks/run.py --datasets syn_100k,syn_500k,syn_2M,syn_6M --tasks 'bins.*,cells.*'
.venv/bin/python benchmarks/report.py                                       # tables + figures
```

Results are appended to `results/<name>.jsonl`, one line per run, holding the phases, memory and the reason for any non-ok status. `report.py` writes the following to `results/report/`:

- `tasks.csv`, `phases.csv` and `outcomes.csv`;
- `summary.md`;
- figures (PDF and PNG): key metrics per tool for Mouse Brain, scaling with dataset size, and task coverage.

Per-run logs and, with `snapshots: true`, a picture of each tool's last rendered view are kept in `data/benchmark_work/runs/` for visual checks.

### On a Slurm compute node

The full Mouse Brain run plus the sweep takes hours, so run it as a batch job on a whole node:

```bash
sbatch --export=ALL,DATASETS=syn_100k,REPEATS=1,NAME=smoke_syn100k benchmarks/slurm/benchmark.sbatch   # quick check
sbatch --export=ALL,DATASETS=mouse_brain,NAME=mouse_brain benchmarks/slurm/benchmark.sbatch
sbatch --export=ALL,DATASETS=syn_100k,syn_500k,syn_2M,syn_6M,NAME=sweep benchmarks/slurm/benchmark.sbatch
```

- `--resume` skips runs already in `results/<NAME>.jsonl`, so resubmitting after a time limit continues where the job stopped.
- Logs go to `benchmarks/slurm/logs/`.
- Submit from the repository root.

## Layout

```
config.yaml          datasets, tools, interpreters, repeats, limits
run.py               orchestrator (one process per run)        report.py   tables + figures
harness/             measure.py (process-tree memory, timeout)  worker.py (phases)
                     browser.py (Playwright settle, frames)     viewports.py, datasets.py, synthetic.py
adapters/            tool_hdviewer.py, tool_hdviewer_browser.py (+ hdviewer_server.py)
                     tool_spatialdata.py, tool_napari.py, tool_vitessce.py (+ vitessce_app/)
                     loupe/ (PROTOCOL.md, loupe_monitor.py, results_template.csv)
envs/                pip requirements of the benchmark environments
```

The earlier viewport-query micro-benchmark is in the repository history (`git show b14081f:benchmarks/viewport_queries.py`).
