# Cell segmentation plugins

Each segmentation method is a self-contained plugin in `plugins/<name>/`. The viewer
lists every plugin it finds on the loading page. Users can select several plugins; each
one runs on its own and becomes a separate cells layer, so the methods can be compared
side by side. The Comparison panel can also match the cells of two methods against each
other.

Plugins run in isolation, normally in a container, so each one keeps its own
dependencies (TensorFlow, PyTorch, …) without affecting the viewer or the other plugins.
Each plugin names a **prebuilt container image** in its manifest. The viewer pulls it
from the registry the first time it is used, so nothing is installed or built at
runtime. The `Dockerfile` next to each plugin is the source of that image. The workflow
in `.github/workflows/segmentation-images.yml` rebuilds and publishes it whenever
`segmentation/` changes on `main`.

| Plugin      | What it does                                            | Needs        |
|-------------|---------------------------------------------------------|--------------|
| `stardist`  | Pretrained H&E nucleus model (`2D_versatile_he`)        | TensorFlow   |
| `cellpose`  | Cellpose `cyto3` (whole cells) or `nuclei`, tiled       | PyTorch      |
| `watershed` | Hematoxylin threshold + watershed. A baseline with no model, and the smallest example to copy | scikit-image |

## The contract

A plugin is a program that:

1. reads `/data/input/image.tif`: the tissue image as `(H, W, 3)` RGB, row 0 at the top,
   in the original resolution of the dataset;
2. reads `/data/input/params.json`:
   ```json
   {"params": {"diameter_um": 10, "model": "cyto3"}, "um_per_px": 0.2738}
   ```
   `params` holds the values from the manifest, after the user's changes and validation.
   `um_per_px` is the physical size of a pixel, or `null` if the dataset does not provide it;
3. writes `/data/output/cells.geojson`: a `FeatureCollection` with one `Polygon` per cell,
   in **pixel coordinates of the input image** (`x` = column, `y` = row). Feature
   `properties` are optional, e.g. `{"score": 0.93}`;
4. optionally prints `##progress <fraction 0-1> <message>` lines, which the loading page
   shows as a progress bar;
5. exits with a non-zero code on failure. The last lines of its output are shown to the
   user and printed in the server console, and the full log is kept.

The input folder is mounted read-only. Write only to the output folder.

The viewer converts the polygons to figure coordinates, assigns each Visium HD bin to the
cell it falls in, and adds the cells as a layer.

## Adding a method

1. Copy `plugins/watershed/` to `plugins/<your_method>/`. The folder name is the plugin id.
2. Edit **`manifest.yaml`**:
   ```yaml
   name: my_method                 # must match the folder name
   title: My method                # shown in the UI
   description: One sentence on what it segments and how.
   version: 0.1.0                  # bump it when results change: cached results are per version
   image: ghcr.io/<owner>/visium-seg-my_method:0.1.0   # prebuilt image, pulled on first use
   url: https://...                # optional
   citation: '...'                 # optional; quote it if it contains ": "
   params:                         # turned into inputs on the loading page
     - name: diameter_um
       label: Expected diameter (µm)
       type: float                 # float | int | bool | choice
       default: 10                 # null = "use the method's own default"
       min: 2
       max: 60
       help: Shown as a tooltip.
     - name: model
       type: choice
       choices: [a, b]
       default: a
   ```
3. Write **`run.py`** using the small SDK in `sdk/vhd_seg.py`:
   ```python
   from vhd_seg import start, iter_blocks, count_blocks, labels_to_polygons

   job = start()                          # parses params.json and --input/--output
   image = job.read_image()
   diameter_px = job.um_to_px(job.params["diameter_um"])

   polygons = []
   for core, (y0, y1, x0, x1) in iter_blocks(image.shape, 2048, overlap=64):
       labels = my_model(image[y0:y1, x0:x1])            # integer mask, 0 = background
       polygons += labels_to_polygons(labels, offset=(y0, x0), core=core)
       job.report(..., "Segmenting")

   job.write_polygons(polygons)
   ```
   Whole-slide images are large (often >20 000 px a side). `iter_blocks` tiles the image
   with overlap. `labels_to_polygons(..., core=core)` keeps each cell from exactly one
   tile, so no cell is duplicated or cut at the tile edges. Methods that already
   produce polygons (like StarDist) can pass them to `job.write_polygons` directly.
4. List the dependencies in **`requirements.txt`** and adapt the **`Dockerfile`**. Bake
   model weights into the image (see `plugins/cellpose/Dockerfile`), so that runs work
   offline and under any user id.

   Instead of a Dockerfile, `image` can point at any existing image that already contains
   your `run.py`, the SDK and their dependencies, with an entrypoint that runs `run.py`.
   The viewer only pulls and runs it.
5. Test it without the viewer:
   ```bash
   cd segmentation
   docker build -t ghcr.io/<owner>/visium-seg-my_method:0.1.0 -f plugins/my_method/Dockerfile .
   mkdir -p /tmp/seg/in /tmp/seg/out
   cp some_image.tif /tmp/seg/in/image.tif
   echo '{"params": {"diameter_um": 10, "model": "a"}, "um_per_px": 0.27}' > /tmp/seg/in/params.json
   docker run --rm -v /tmp/seg/in:/data/input:ro -v /tmp/seg/out:/data/output ghcr.io/<owner>/visium-seg-my_method:0.1.0
   ```
   You can also run it without Docker:
   `PYTHONPATH=sdk python plugins/my_method/run.py --input /tmp/seg/in --output /tmp/seg/out`.

6. Publish the image: merging to `main` runs the workflow, which pushes it to the tag in
   `manifest.yaml`. Packages on ghcr.io are private by default: make each one **public**
   once (GitHub → Packages → package settings → Change visibility), so the viewer can
   pull it without logging in.

Restart the viewer and the new method appears on the loading page. No viewer code
needs to change.

Always use a registry namespace you control (here `ghcr.io/luisalonsoesteban/`).
The viewer pulls and runs whatever image the manifest names, so a name under someone
else's account would run their code.

## Runners

`config.yaml` (or the `VHD_SEG_RUNNER` environment variable) chooses how plugins run:

- **`docker`**: `docker run` with the input and output folders mounted. The image is
  pulled from its registry the first time it is needed. If it is not published (a
  plugin under development), it is built from the plugin's `Dockerfile` instead. Set
  `gpu: true` to pass `--gpus all`.
- **`apptainer`**: for HPC clusters without Docker. Runs `docker://<image>`, which
  Apptainer pulls and caches as a SIF on first use. To use a local `.sif` instead
  (e.g. on nodes without internet access), put it at `images/<name>.sif` or set its
  path in `apptainer.images.<name>`:
  `apptainer pull images/cellpose.sif docker://ghcr.io/luisalonsoesteban/visium-seg-cellpose:0.1.0`.
- **`local`**: runs `run.py` with a Python interpreter on this machine. By default
  that is the viewer's own interpreter. Point a plugin at its own environment with
  `local.python.<name>: /path/to/env/bin/python`. Useful for development, or where no
  container runtime is available.
- **`auto`** (default): docker if the daemon is reachable, else apptainer, else local.

`max_parallel` (default 2) is how many methods run at the same time; further selected
methods wait in a queue (shown as "Queued" on the loading page) and start as soon as one
finishes. Segmentation is memory hungry, so raise it only with enough RAM for that many
runs on a large slide. `threads_per_run` (default 8) caps the CPU threads of each run through `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, the TensorFlow pool variables and Docker's
`--cpus`. Without a cap, numerical libraries start one thread per core in every run,
which exceeds per-user thread limits on shared servers. Plugins that manage their own
threads can read the budget from `VHD_SEG_THREADS`.

## Results and caching

Runs are stored next to the dataset, in `<dataset dir>/segmentation/<plugin>/<hash>/`,
together with `run.log`. The hash covers the image, the plugin version and the parameters.
Loading the same dataset with the same settings reuses the result instantly. Failed runs
stay in `<hash>.running/` with their log, for debugging. The first line of `run.log` is the
exact command that was run, and the last line says why the run ended. A run the kernel
kills for exceeding a memory limit cannot print an error itself, so the runner writes
`killed, most likely out of memory` there and on the loading page.
