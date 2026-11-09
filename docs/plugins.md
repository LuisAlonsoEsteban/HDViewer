# Plugins: bring your own method

VisiumVisualizer runs its **cell segmentation** and **cell-type deconvolution**
methods as plugins. A plugin is one folder holding a short manifest and a script.
The viewer finds it at start-up and builds its settings form from the manifest:
inputs, defaults, limits and help texts. It runs the method in isolation and
shows the result as a layer or a readout. **No viewer code changes**: if your
method reads an image or a count matrix and writes polygons or proportions, it can
be a plugin.

This page is the starting point for method developers. Each plugin kind has a full
reference next to its code:

- Segmentation: [`hd_viewer/startup/data/image/segmentation/api/README.md`](../hd_viewer/startup/data/image/segmentation/api/README.md)
- Deconvolution: [`hd_viewer/analysis/deconvolution/api/README.md`](../hd_viewer/analysis/deconvolution/api/README.md)

## The two kinds of plugin

| | Segmentation | Deconvolution |
|---|---|---|
| Where it appears | Loading screen, under *Cell segmentation* | Each layer's **Deconv.** tab |
| What it receives | The tissue image (`image.tif`, RGB) and its µm per pixel | The layer's counts (`target.h5ad`) and the single-cell reference (`reference.h5ad`) |
| What it returns | One polygon per cell (`cells.geojson`, image pixels) | Cell-type fractions per observation (`fractions.json`), optional abundances |
| What the viewer does with it | Assigns the bins to the cells and adds a **Cells · your method** layer | Stores the run with the layer, colours the map by any cell type, and offers it to Groups and Metadata comparison |
| Plugin folder | `hd_viewer/startup/data/image/segmentation/api/plugins/` | `hd_viewer/analysis/deconvolution/api/plugins/` |
| Helper module (SDK) | `sdk/vhd_seg.py` | `sdk/vhd_deconv.py` |
| Simplest example to copy | `plugins/watershed/` (no model, ~60 lines) | `plugins/nnls/` |

Every plugin folder contains:

```text
plugins/my_method/
├── manifest.yaml      # name, title, version, container image, parameters
├── run.py             # reads the inputs, runs the method, writes the result
├── requirements.txt   # the method's own dependencies
└── Dockerfile         # builds the image named in the manifest (optional, see below)
```

## How plugins run

The viewer never imports a plugin. For every run it writes the inputs to a folder,
starts the plugin and reads its result files:

1. The input folder is mounted read-only at `/data/input`, and an empty output
   folder at `/data/output`. Outside a container, `run.py` gets `--input` and
   `--output` arguments instead; the SDK handles both cases.
2. The plugin prints `##progress <fraction 0-1> <message>` lines, which the
   interface shows as a progress bar. The SDK's `job.report()` does this.
3. A non-zero exit code is a failure: the last lines of the output are shown to the
   user, and the full log is kept with the run.

A **runner** starts the plugin. It is configured per service in `config.yaml`, or
with `VHD_SEG_RUNNER` / `VHD_DECONV_RUNNER`:

- `docker`: one container per run. The manifest's image is pulled on first use, or
  built from the plugin's `Dockerfile` if it is not published yet.
- `apptainer`: the same images on HPC clusters without Docker.
- `local`: `run.py` with a Python environment on the same machine. This is the
  quickest way to develop a plugin.
- `auto` (default): Docker if available, else Apptainer, else local.

Each method keeps its own dependencies (TensorFlow, PyTorch, R…) in its own image or
environment, so plugins never conflict with the viewer or with each other. Results
are cached per input, plugin **version** and parameters: rerunning the same
settings is instant, and bumping `version` invalidates old results.

## Quick start: a segmentation method

1. **Copy the template.** From the repository root:

   ```bash
   cd hd_viewer/startup/data/image/segmentation/api
   cp -r plugins/watershed plugins/my_method
   ```

2. **Edit `plugins/my_method/manifest.yaml`.** `name` must match the folder:

   ```yaml
   name: my_method
   title: My method
   description: One sentence on what it segments and how.
   version: 0.1.0
   image: ghcr.io/<your-account>/visium-seg-my_method:0.1.0
   url: https://github.com/<you>/my-method
   citation: 'You et al. (2026) My method. Journal.'
   params:
     - name: diameter_um
       label: Expected diameter (µm)
       type: float
       default: 10
       min: 2
       max: 60
       help: Typical cell diameter.
   ```

3. **Write `run.py`** with the SDK. Whole-slide images are large, so process them in
   overlapping blocks. `labels_to_polygons(..., core=core)` keeps each cell exactly
   once:

   ```python
   from vhd_seg import count_blocks, iter_blocks, labels_to_polygons, start

   job = start()                               # reads params.json; parses --input/--output
   image = job.read_image()                    # (H, W, 3) RGB, row 0 at the top
   diameter = job.um_to_px(job.params["diameter_um"])

   total, polygons = count_blocks(image.shape, 2048), []
   for i, (core, (y0, y1, x0, x1)) in enumerate(iter_blocks(image.shape, 2048, overlap=64)):
       labels = my_model(image[y0:y1, x0:x1], diameter)   # integer mask, 0 = background
       polygons += labels_to_polygons(labels, offset=(y0, x0), core=core)
       job.report(i / total, f"Segmenting tile {i + 1}/{total}")

   job.write_polygons(polygons)                # cells.geojson in image pixel coordinates
   ```

   Methods that already return polygons (as StarDist does) can pass them to
   `job.write_polygons` directly.

4. **Test it without the viewer**, on any RGB TIFF:

   ```bash
   mkdir -p /tmp/seg/in /tmp/seg/out
   cp some_image.tif /tmp/seg/in/image.tif
   echo '{"params": {"diameter_um": 10}, "um_per_px": 0.27}' > /tmp/seg/in/params.json
   PYTHONPATH=sdk python plugins/my_method/run.py --input /tmp/seg/in --output /tmp/seg/out
   ```

5. **Run it in the viewer.** Restart the server (`uvicorn hd_viewer.api.main:app`).
   *My method* now appears on the loading screen, and ticking it adds a
   **Cells · My method** layer next to the other methods, ready to compare. During
   development, set `VHD_SEG_RUNNER=local`, and point `local.python.my_method` in
   `config.yaml` at an environment with your `requirements.txt` installed.

## Quick start: a deconvolution method

1. **Copy the template.** From the repository root:

   ```bash
   cd hd_viewer/analysis/deconvolution/api
   cp -r plugins/nnls plugins/my_deconv
   ```

2. **Edit `plugins/my_deconv/manifest.yaml`** (same keys as above). Parameters of type
   `obs_column` let the user pick a column of the reference, for example a donor
   annotation.

3. **Write `run.py`.** The inputs are AnnData files with untransformed counts in `.X`;
   `job.label_column` is the reference column holding the cell types. This complete,
   if simplistic, example assigns each observation to its most similar cell type:

   ```python
   import numpy as np
   import pandas as pd
   from vhd_deconv import start

   job = start()
   target, reference = job.read_inputs()
   genes = target.var_names.intersection(reference.var_names)
   labels = reference.obs[job.label_column].astype(str).to_numpy()

   # Mean profile of each cell type, as proportions of its counts
   profiles = reference[:, genes].to_df().groupby(labels).mean()
   profiles = profiles.div(profiles.sum(axis=1), axis=0)

   counts = target[:, genes].to_df()
   totals = counts.sum(axis=1)
   assigned = (totals > 0).to_numpy()          # all-zero observations stay unassigned (NaN)
   job.report(0.5, "Scoring cell types")
   similarity = counts[assigned].div(totals[assigned], axis=0) @ profiles.T
   fractions = pd.DataFrame(np.nan, index=target.obs_names, columns=profiles.index)
   fractions.loc[assigned] = np.eye(len(profiles))[similarity.to_numpy().argmax(axis=1)]

   job.write_result(fractions, genes)          # fractions.json + result.json
   ```

   The viewer checks the result before storing it. Rows must match `target.obs_names`
   in order, and columns must match the reference's cell types. Each assigned row must
   be non-negative and sum to one; unassigned rows are entirely NaN. Absolute
   estimates (e.g. cells per spot) go in `abundance=`, with their `abundance_units`.

4. **Test it without the viewer:**

   ```bash
   mkdir -p /tmp/deconv/in /tmp/deconv/out
   cp target.h5ad reference.h5ad /tmp/deconv/in/
   echo '{"params": {}, "label_column": "cell_type"}' > /tmp/deconv/in/params.json
   PYTHONPATH=sdk python plugins/my_deconv/run.py --input /tmp/deconv/in --output /tmp/deconv/out
   ```

5. **Run it in the viewer.** Restart the server and open a layer's **Deconv.** tab:
   your method is in the list, with a form built from its manifest. The reference is
   uploaded once on the loading screen and reused by every layer.

Methods that fit a reference signature, as cell2location and CIBERSORTx do, can
cache it across layers through `reference_params`, `job.reference_cache` and
`job.reference_output`. See the deconvolution reference for details.

## Manifest reference

| Key | Required | Meaning |
|---|---|---|
| `name` | yes | Plugin id; must equal the folder name |
| `title` | yes | Name shown in the interface |
| `image` | yes | Container image of the method (`registry/owner/name:tag`) |
| `version` | | Cache key: bump it whenever results change (default `0`) |
| `description` | | One or two sentences shown above the form |
| `url` | | Documentation link shown with the method |
| `citation` | | Reference shown with the method ("Please cite: …"); quote it in YAML |
| `params` | | The settings form (see below) |
| `runners` | | Runners the method supports (default `[docker, apptainer, local]`) |
| `command`, `input_mount`, `output_mount` | | Override the container command and mount points, e.g. to wrap an existing image |
| `reference_params` | | Deconvolution: parameters that change the fitted reference signature |

Each entry of `params` becomes one input:

| Field | Meaning |
|---|---|
| `name` | Key in `job.params` |
| `type` | `float`, `int`, `bool`, `choice` (with `choices`), or `obs_column` (deconvolution: a reference column) |
| `default` | Initial value; `null` means "the method's own default" |
| `label`, `help` | Text of the input and its description |
| `min`, `max`, `step` | Limits of numeric inputs, enforced by the viewer |
| `required` | The user must provide a value |

Invalid manifests are reported in the server log and skipped, so the other methods
keep working.

## Packaging the method in a container

Containers make a method reproducible, and are how the public instance runs it.

- **Build from the bundled pattern.** Copy the template's `Dockerfile`. It installs
  `requirements.txt`, then copies the SDK and `run.py` into `/app`. Bake model
  weights into the image (see [cellpose's Dockerfile](../hd_viewer/startup/data/image/segmentation/api/plugins/cellpose/Dockerfile)), so
  runs work offline and under any user id. Build from the service's `api/` folder:

  ```bash
  docker build -t ghcr.io/<you>/visium-seg-my_method:0.1.0 -f plugins/my_method/Dockerfile .
  docker run --rm -v /tmp/seg/in:/data/input:ro -v /tmp/seg/out:/data/output ghcr.io/<you>/visium-seg-my_method:0.1.0
  ```

- **Or wrap an existing image.** If your method already ships a container, name it
  in `image`, and use `command` (and `input_mount` / `output_mount` if needed) to call
  a small `run.py` that reads and writes the files described above.

- **Publish it** in a registry you control, such as `ghcr.io/<your-account>/…`, and
  make the package public so the viewer can pull it without logging in. The
  repository's GitHub Actions workflows (`.github/workflows/*-images.yml`) build and
  publish every plugin folder that has a `Dockerfile`.

## Sharing your method

There are two ways to share a method:

- **Use it in your own copy.** Keep the plugin folder in your fork or local clone.
  Anyone who runs VisiumVisualizer from it gets your method. This is enough for a lab,
  a paper's reproducibility package, or a private server.
- **Add it to VisiumVisualizer and <https://hdviewer.unav.edu>.** Open a pull request
  that adds the plugin folder; the [contribution guide](../CONTRIBUTING.md) has the
  checklist. Once it is merged and its image is published, the method is available to
  everyone who uses the public instance, with your description, link and citation.

The public server only runs plugins merged into this repository. Visitors cannot
upload code, so every method there has been reviewed. If you are unsure whether your
method fits, open a [*New method plugin* issue](https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues/new?template=method-plugin.yml)
first.

### Before you share: checklist

- [ ] `name` matches the folder; `title`, `description`, `url` and `citation` are filled in.
- [ ] Every parameter has a `label`, sensible `default`, `min`/`max` where relevant, and `help`.
- [ ] `run.py` reports progress, fails with a clear message on bad input, and writes only to the output folder.
- [ ] The image builds from the `Dockerfile` (or `image` names a public image) and runs without network access once pulled.
- [ ] Results are deterministic for fixed parameters (set seeds); bump `version` whenever they change.
- [ ] The method's licence allows redistribution in a container, and is stated in its folder or image.
- [ ] Tested standalone and in the viewer on an example dataset (`python -m unittest discover -s tests` still passes).

## Security

Running a plugin means running its code and image. Only add plugins from sources you
trust. Always name images under a registry account you control: the viewer pulls
and runs exactly the image a manifest names. Never put credentials in a manifest or
in parameters. Methods that need a licence key, like CIBERSORTx, read it from the
server's environment variables, and the runner keeps it out of logs and caches.
