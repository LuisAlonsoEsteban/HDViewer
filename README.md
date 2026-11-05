# VisiumVisualizer

An interactive web platform for exploring **10x Genomics Visium HD** spatial
transcriptomics data on top of the original tissue image: bins, pseudospots and
segmented cells; gene expression, clustering, differential expression, regions of interest, cell-type deconvolution and sample comparison.

**Use it online, no installation or login required:** <https://hdviewer.unav.edu>

The public server is free and open to all users. Uploaded data is private to your browser session and is deleted from the server when the session ends. For large datasets, private data you prefer not to upload, or offline work, you can run the same application on your own computer; the rest of this document explains how.

- [Features](#features)
- [Run it locally](#run-it-locally)
- [Example datasets](#example-datasets)
- [Input files](#input-files)
- [Optional: segmentation and deconvolution](#optional-segmentation-and-deconvolution)
- [Configuration](#configuration)
- [Hosting a public instance](#hosting-a-public-instance)
- [Documentation](#documentation) · [Data privacy](#data-privacy) · [Citation](#citation) · [License](#license) · [Contact](#contact)

## Features

- **Multi-resolution viewing** of Visium HD bins (2, 8, 16 µm and aggregates),
  Visium v2-like 55 µm pseudospots and segmented cells, over a tiled
  full-resolution H&E image.
- **Expression and metadata maps** with gene/metadata autocomplete, genes shown
  as raw, normalized or log-normalized counts, and expression-based queries.
- **Clustering** (Leiden, k-means, DBSCAN, HDBSCAN, Gaussian mixture,
  agglomerative) with marker genes and UMAP/PCA embeddings.
- **Groups and regions of interest** drawn on the tissue, by radius around seeds,
  or by barcode, which can be compared by **differential expression** (Wilcoxon
  rank-sum test, Benjamini–Hochberg correction).
- **Cell-type deconvolution** against a single-cell reference with NNLS,
  RCTD (spacexr), cell2location, MuSiC or CIBERSORTx.
- **Sample comparison** across layers or resolutions.
- Export of layer data as AnnData (`.h5ad`) and of the tissue map as PNG.

## Run it locally

### Requirements

- Linux, macOS or Windows (WSL recommended on Windows).
- Python **3.11** or newer.
- Memory, as a rough guide: 16 GB RAM for 8/16 µm bins; 32–64 GB for 2 µm bins of
  a full slide or cell segmentation.
- A recent desktop browser: Chrome, Firefox, Edge or Safari.
- Optional: Docker or Apptainer/Singularity for cell segmentation and most
  deconvolution methods (see [below](#optional-segmentation-and-deconvolution)).

### Install

```bash
git clone https://github.com/LuisAlonsoEsteban/VisiumVisualizer.git
cd VisiumVisualizer
python3 -m venv .venv
source .venv/bin/activate          # Windows (PowerShell): .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Using conda instead:

```bash
conda create -n visiumvisualizer python=3.11
conda activate visiumvisualizer
pip install -r requirements.txt
```

### Launch

From the repository root:

```bash
bokeh serve hd_viewer --show
```

Your browser opens <http://localhost:5006/hd_viewer>. The landing page links
to the tutorial; **Launch** opens the loading screen, where you either pick an
example dataset or upload your own files.

- Uploads stream to the server in 4 MB pieces, so multi-gigabyte Visium HD files
  need no extra flags. Each file may be up to 20 GiB (`VHD_UPLOAD_MAX_BYTES`), and
  the `tmp/` folder in the repository needs enough free disk space to hold it.
- Image tiles are served by a small helper server on port **4000**. If that port is
  taken, choose another with `VHD_TILE_PORT=4100 bokeh serve hd_viewer ...`.
- Use `--port 5007` to change the viewer's port.

### Running on a remote machine (HPC, workstation)

Start the viewer on the remote machine without `--show`, then forward both ports
from your laptop and open <http://localhost:5006/hd_viewer> locally:

```bash
ssh -L 5006:localhost:5006 -L 4000:localhost:4000 user@remote-host
```

## Example datasets

The loading screen offers every example dataset whose files are present on the
machine. Datasets are looked up under `data/` in the repository (or the folder in
`VHD_DATA_DIR`), with this layout:

```text
data/
├── Visium_HD_Mouse_Brain/
│   ├── Visium_HD_Mouse_Brain_tissue_image.tif
│   └── binned_outputs/square_{002,008,016}um/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/tissue_positions.parquet
└── CRC_P1/
    ├── VisiumHD/
    │   ├── Visium_HD_Human_Colon_Cancer_P1_tissue_image.btf
    │   └── binned_outputs/square_008um/
    │       ├── filtered_feature_bc_matrix.h5
    │       └── spatial/tissue_positions.parquet
    └── Chromium/                                   # single-cell reference, optional
        ├── HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5
        └── SingleCell_MetaData.csv
```

Both datasets are public:

- **Human colorectal cancer, patient P1** (Oliveira et al., *Nature Genetics*, 2025,
  [doi:10.1038/s41588-025-02193-3](https://doi.org/10.1038/s41588-025-02193-3)),
  with its matched Chromium Flex reference:
  <https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc>.
- **Mouse brain (Visium HD, FFPE)**: search for "Visium HD mouse brain" in the
  [10x Genomics datasets portal](https://www.10xgenomics.com/datasets).

Download the "binned outputs" archive and the full-resolution tissue image for
each, extract them into the layout above, and restart the viewer. To add your own
examples, edit [`hd_viewer/startup/datasets.py`](hd_viewer/startup/datasets.py).

## Input files

Upload the standard Space Ranger outputs for **one bin size** (2, 8 or 16 µm):

| File | Where to find it in the Space Ranger output | Required |
|---|---|---|
| Expression matrix (`.h5`) | `binned_outputs/square_XXXum/filtered_feature_bc_matrix.h5` | yes |
| Bin positions (`.parquet`) | `binned_outputs/square_XXXum/spatial/tissue_positions.parquet` | yes |
| Tissue image (`.btf`, `.tif`, `.tiff`, `.png`) | the full-resolution image given to Space Ranger | yes |
| Single-cell reference (`.h5ad` or 10x `.h5`) | your annotated scRNA-seq data, for deconvolution | no |
| Reference annotations (`.csv`, `.parquet`) | barcode in the first column, then cell-type labels | no |

The [tutorial](hd_viewer/static/tutorial.html) explains each option on the loading
screen and how to interpret every result.

## Optional: segmentation and deconvolution

Cell segmentation (Cellpose, StarDist, watershed) and the deconvolution methods
run as plugins in containers, so their dependencies do not have to be installed
alongside the viewer. The runner is picked automatically: Docker if its daemon is
reachable, otherwise Apptainer/Singularity, otherwise the local Python environment.
Images are pulled from `ghcr.io` on first use.

- Segmentation: [`hd_viewer/startup/data/image/segmentation/api/config.yaml`](hd_viewer/startup/data/image/segmentation/api/config.yaml)
- Deconvolution, including RCTD/MuSiC R setups and CIBERSORTx credentials
  (`VHD_CIBERSORTX_USERNAME`, `VHD_CIBERSORTX_TOKEN`):
  [deconvolution plugin guide](hd_viewer/board/operations/layer/deconvolution/api/README.md)

NNLS deconvolution and every other feature work without containers.

Each run is limited in time (`timeout_minutes`: 120 for segmentation, 720 for
deconvolution, as RCTD on a full slide takes about 8 hours) and, with Docker, in
memory (`memory`: 32g and 64g) and processes (`pids_limit`). Set these in the two
`config.yaml` files above; an empty value disables a limit. Apptainer applies
`apptainer.memory` only on hosts with cgroups v2 delegation. Runs belong to the
browser session that started them: closing the tab cancels them, unless another
session is waiting for the same result. If a run fails, users see a short message
with a run id; its log is kept as `<run id>/run.log` under the dataset's
`segmentation/` or `deconvolution/` folder.

## Configuration

Defaults live in [`hd_viewer/settings/__init__.py`](hd_viewer/settings/__init__.py).
The following environment variables override them without editing code:

| Variable | Default | Purpose |
|---|---|---|
| `VHD_DATA_DIR` | `<repo>/data` | Folder holding the example datasets |
| `VHD_TILE_PORT` | `4000` | Port of the image tile server |
| `VHD_TILE_HOST` | `0.0.0.0` | Address the tile server binds to |
| `VHD_TILE_URL` | `http://localhost:<port>` | Tile URL as seen by the browser (see below) |
| `VHD_TILE_FORMAT`, `VHD_TILE_QUALITY` | `webp`, `80` | Tile encoding (`webp` or lossless `png`); changing it re-tiles |
| `VHD_UPLOAD_MAX_BYTES` | 20 GiB | Largest file accepted per upload field |
| `VHD_SEG_RUNNER`, `VHD_DECONV_RUNNER` | `auto` | Force `docker`, `apptainer` or `local` plugins |
| `VHD_SEG_THREADS`, `VHD_DECONV_THREADS` | `8` | CPU threads per segmentation/deconvolution run |
| `VHD_DECONV_GPU` | off | Pass GPUs to deconvolution containers |

Uploaded files and per-session results are written to `tmp/<session>/` in the
repository and removed when the browser session closes.

## Hosting a public instance

This is how <https://hdviewer.unav.edu> is deployed; you only need it to host your
own shared instance. Run a single Bokeh process (the tile registry is per process)
behind an HTTPS reverse proxy that also routes the tile server:

```bash
export VHD_TILE_HOST=127.0.0.1 VHD_TILE_URL=/tiles
export BOKEH_SECRET_KEY=$(python -c "import secrets; print(secrets.token_urlsafe(32))")
export BOKEH_SIGN_SESSIONS=yes
bokeh serve hd_viewer --port 5006 --address 127.0.0.1 \
    --allow-websocket-origin=hdviewer.unav.edu --use-xheaders
```

```nginx
server {
    listen 443 ssl;
    server_name hdviewer.unav.edu;
    location /tiles/ {
        proxy_pass http://127.0.0.1:4000/;
    }
    location / {
        proxy_pass http://127.0.0.1:5006;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }
}
```

## Documentation

- **Tutorial:** [`hd_viewer/static/tutorial.html`](hd_viewer/static/tutorial.html),
  also at `/hd_viewer/static/tutorial.html` on a running instance and behind every
  in-app help button.
- **Biological use case:** a walkthrough reproducing Figure 4a–b of Oliveira et al.
  (2025) on the P1 colorectal cancer sample:
  [`hd_viewer/static/biological-use-case.html`](hd_viewer/static/biological-use-case.html).
- **Plugin API** for adding segmentation or deconvolution methods:
  [deconvolution plugins](hd_viewer/board/operations/layer/deconvolution/api/README.md).

## Data privacy

Each browser session gets its own private folder on the server. Uploaded files,
intermediate results and exports are only reachable from that session, and the
folder is deleted when the session ends (closing the tab, or after inactivity).
Exported `.h5ad` files get a random, unguessable name and are removed from the
server ten minutes after they are prepared. No account, e-mail address, analytics or tracking cookies are used.

## Running the tests

```bash
python -m unittest discover -s tests
```

## Citation

If you use VisiumVisualizer, please cite the software (see
[`CITATION.cff`](CITATION.cff); GitHub's "Cite this repository" button exports it)
and the methods and datasets you used in your analysis.

## License

VisiumVisualizer is released under the [MIT License](LICENSE). Example datasets and
third-party methods run as plugins keep their own licenses.

## Contact

Luis Alonso, Computational Biology group, Tecnun, University of Navarra —
<laesteban@unav.es>. Please report bugs and feature requests on the
[issue tracker](https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues).

---

## For developers

Presentation files:

- `hd_viewer/templates/index.html`: landing page and workspace shell.
- `hd_viewer/static/css/app.css`: responsive landing/setup/workspace layout.
- `hd_viewer/static/css/widgets.css`: shared styles inside Bokeh shadow roots.
- `hd_viewer/presentation.py`: project metadata and per-session style integration.
- `hd_viewer/static/assets/README.md`: image/logo sources and attribution.

`CITATION.cff` and its downloadable copy at `hd_viewer/static/citation.cff` must be
kept identical.

The biological walkthrough lives in `hd_viewer/static/biological-use-case.html`.
Give each new section a unique `id` and add its link to the table of contents in
that file. For screenshots, create `hd_viewer/static/use-cases/` and reference
images with relative paths such as `use-cases/result.png`. Shared styling lives in
`hd_viewer/static/tutorial.css`.

Performance benchmarks are described in [`benchmarks/README.md`](benchmarks/README.md).
