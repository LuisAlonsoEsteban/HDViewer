# Loupe Browser: manual measurement protocol

Loupe Browser (10x Genomics) is a closed-source desktop application for Windows and macOS. It has no command-line or scripting interface, so it cannot run in the automated suite. Its tasks are measured by hand with `loupe_monitor.py`, which:

- records the time between two key presses;
- samples the memory of every Loupe Browser process at 10 Hz.

The rows it writes use the same task and phase names as the automated tools. `benchmarks/report.py --loupe <csv>` merges them into the same tables and figures.

## Setup

1. Use a machine comparable to the benchmark node, and record its CPU, RAM, OS and GPU in the `notes` of the first row.
2. Install the current Loupe Browser and record its version.
3. Install the monitor: `pip install psutil`. Copy `loupe_monitor.py` to the machine.
4. Get the input. Loupe opens the `.cloupe` file that Space Ranger writes: `data/Visium_HD_Mouse_Brain/binned_outputs/square_008um/cloupe.cloupe`. Copy it to a local SSD.
5. Close other applications. Start screen recording if you want frame-accurate times (see *Accuracy*).

## Procedure

Run `python loupe_monitor.py --dataset mouse_brain --out loupe_results.csv`. For each row below:

- type the task, phase and bin size at the prompt;
- press Enter as you click;
- press Enter again when the view has finished drawing: no spinner, and every tile and spot is sharp.

Repeat each row **5 times**. Restart Loupe Browser before each repeat of `e2e.first_view`.

| Task | Phase | Bin | Action in Loupe Browser |
|---|---|---|---|
| e2e.first_view | load_to_first_view | 8 | File › Open `cloupe.cloupe` → spatial view shown with image and spots |
| image.render | view:overview | | Hide spots (Spatial › opacity 0), click *fit to window* |
| image.render | view:zoom_250um | | Zoom until the scale bar reads ~250 µm across the view, centred on tissue |
| image.render | view:zoom_50um | | Zoom to ~50 µm across |
| image.render | pan | | Ten drags of half the view width (right ×2, down ×2, left ×3, up ×2, right ×1); press Enter after the last one has drawn |
| bins.render | colour_by_gene | 8 | Show spots, then type `Ttr` in *Feature* search and select it |
| bins.render | view:overview | 8 | Fit to window |
| bins.render | view:zoom_1000um | 8 | Zoom to ~1 mm across |
| bins.render | view:zoom_250um | 8 | Zoom to ~250 µm across |
| bins.render | view:zoom_50um | 8 | Zoom to ~50 µm across |
| bins.render | pan | 8 | The ten drags above at ~250 µm |
| bins.render | gene_switch | 8 | Replace `Ttr` with `Mbp` |
| bins.de | total | 8 | Locally distinguishing features between two clusters (if offered) |

The following rows are recorded as `n.a.` (type them with the note `n.a.: <reason>`, then press Enter twice):

- `image.pyramid`, `bins.load` and `cells.*`: precomputed by Space Ranger, or not importable. Loupe cannot import external cell outlines such as `cells.geojson`.
- `bins.query` and `bins.aggregate`: no user-facing equivalent.
- `bins.cluster`: Loupe shows Space Ranger's precomputed clusters.
- 2 and 16 µm bins: Loupe shows the bin sizes stored in the `.cloupe` file. Record each one Loupe offers by switching the bin size in the spatial view (bin column `2` or `16`).

## Accuracy

Key presses add the operator's reaction time, roughly 0.2–0.3 s on each press, with less variance than its size. Report Loupe's figures as manual measurements.

For frame-accurate times:

1. Screen-record at 60 fps.
2. Read the frame of the click and the frame where drawing ends.
3. Overwrite `seconds` in the CSV.

Memory is the summed resident memory of all Loupe processes, as for the browser-based tools.

## Merge

```bash
.venv/bin/python benchmarks/report.py --loupe loupe_results.csv
```
