# Deconvolution plugins

The Deconvolution tab discovers `plugins/<name>/manifest.yaml`. Choose one method,
configure it, then run it. Saved runs remain selectable and export with the layer.
Upload the annotated single-cell reference once on the loading screen. Both inputs
must have unique gene identifiers and nonnegative linear counts in `.X`; matching
names identify shared genes. Methods do not use the current viewport or bin display
aggregation. All-zero observations remain unassigned.

| Method | Additional input | Readouts |
| --- | --- | --- |
| NNLS (default) | Cell-type labels | Expression mixture fractions |
| cell2location | Explicit expected cells/location; optional reference batch | Posterior mean expected cells and fractions normalized from those means |
| MuSiC | Donor/sample annotation with at least two subjects, distinct from cell types | Weighted proportions from `music_prop` |
| spacexr (RCTD) | Raw integer counts; at least two reference cell types | Full, doublet, or multi-mode fractions after platform-effect normalization |
| CIBERSORTx | Container runtime and server credentials | Relative fractions; signature generated from the reference |

Fractions sum to one for assigned rows. These methods have different statistical
interpretations; normalized cell2location means are not a posterior mean of a
ratio. At coarse bin levels fractions average and expected cell abundance sums.
MuSiC and CIBERSORTx are bulk-mixture methods applied to each spatial observation.

## Execution and setup

`config.yaml` controls deconvolution independently of segmentation. Application
settings supply defaults. `VHD_DECONV_RUNNER` overrides `runner`, which may be
`auto`, `docker`, `apptainer`, or `local`. Automatic selection tries reachable
Docker, Apptainer/Singularity, then local Python. Configure one runner for this
service; a plugin cannot silently substitute another method or runner.

`threads_per_run` defaults to 8; `max_parallel` defaults to 1 across layers in a
viewer process. Numerical thread budgets reach plugins via BLAS/OpenMP variables
and `VHD_DECONV_THREADS`. `gpu: true` exposes container GPUs and permits
cell2location to use CUDA. Local cell2location also respects this setting.

NNLS works with the viewer's Python dependencies. Other methods keep their
requirements in their own plugin directory. For local execution, create an
isolated environment, install that method's `requirements.txt`, and add its Python
path to `local.python`:

```yaml
runner: local
local:
  python:
    cell2location: /path/to/cell2location/bin/python
    music: /path/to/music/bin/python
    spacexr: /path/to/spacexr/bin/python
```

MuSiC additionally needs `Rscript` on PATH (or `VHD_MUSIC_RSCRIPT` pointing to it),
the official `xuranw/MuSiC` R package, `SingleCellExperiment`, `Matrix`, and
`jsonlite`. The bundled Dockerfile installs its R dependencies. Sparse reference
counts cross the Python/R boundary as Matrix Market files; only bounded target
batches become dense in R.

NNLS, cell2location, MuSiC, and spacexr pull versioned images from GHCR, using the same
registry namespace as segmentation:

- `ghcr.io/luisalonsoesteban/visium-deconv-nnls:0.1.0`
- `ghcr.io/luisalonsoesteban/visium-deconv-cell2location:0.1.0`
- `ghcr.io/luisalonsoesteban/visium-deconv-music:0.1.0`
- `ghcr.io/luisalonsoesteban/visium-deconv-spacexr:0.1.0`

Docker uses an existing local image or pulls the manifest image. As with segmentation,
a bundled Dockerfile remains a fallback if pulling fails. Apptainer uses a configured
`apptainer.images` entry or local `images/<plugin-name>.sif` when provided; otherwise
it pulls `docker://<manifest image>` and caches the SIF. Apptainer does not build images.

The `Deconvolution plugin images` GitHub Actions workflow builds and publishes these
four images when their sources change on `main`, or when run manually. The workflow
must complete successfully and GHCR packages must be public or accessible to the
execution host before remote pulls can succeed. Changing manifests alone does not
publish images. Bump both the manifest version and image tag when changing methods.
CIBERSORTx retains its official container rather than republishing it.

To build an image manually, run from this API directory:

```bash
docker build -t ghcr.io/luisalonsoesteban/visium-deconv-music:0.1.0 -f plugins/music/Dockerfile .
docker build -t ghcr.io/luisalonsoesteban/visium-deconv-spacexr:0.1.0 -f plugins/spacexr/Dockerfile .
```

### spacexr (RCTD)

The adapter follows the upstream [RCTD vignette](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html):
construct `Reference` and `SpatialRNA`, then call `create.RCTD` and `run.RCTD`.
Local execution needs `Rscript` on PATH (or `VHD_SPACEXR_RSCRIPT`), `Matrix`,
`jsonlite`, and the **GitHub** `dmcable/spacexr` package with this API. The
Bioconductor package has a different API. Install into the R library visible to
that Rscript, for example:

```r
install.packages(c("remotes", "jsonlite"), repos="https://cloud.r-project.org")
remotes::install_github("dmcable/spacexr", build_vignettes=FALSE, upgrade="never")
```

No donor/sample column is required. `.X` must contain untransformed integer
counts; fractional or normalized expression fails validation. The default
reference filter requires 100 UMIs per cell and 25 cells per type after filtering.
RCTD also selects marker genes and filters target depth; sparse HD bins may
need aggregation or adjusted cutoffs. `max_cores` cannot exceed the runner's
`threads_per_run`. A seed controls reference downsampling.

Choose `doublet` (default, at most two types), `full` (all types), or `multi`
(up to `MAX_MULTI_TYPES`). Full and multi weights are normalized per observation.
In doublet mode, singlets receive a unit fraction for their first type;
`doublet_certain` uses the two fitted weights. Rejected and `doublet_uncertain`
observations remain wholly unassigned: an uncertain second identity does not
define confident mixture fractions. Filtered observations retain their original
rows as unassigned. Native confidence information remains in `rctd_results.rds`
and `r_inputs/rctd.json`; `cell_type_labels.json` maps R's safe identifiers back
to the original labels, including slashes and numeric strings.

Counts cross into R as sparse Matrix Market files. RCTD fits the whole layer
because its platform-effect and noise estimates depend on the dataset; arbitrary
independent batches would change the model. Upstream fitting can create dense
matrices, so large layers need sufficient RAM. No fitted model is reused between
layers, since its normalization depends on the target counts.

When available, `target.obsm['spatial']` (an observations-by-2 matrix) is retained.
Otherwise `SpatialRNA` uses its documented placeholder-coordinate option, with
a progress message. RCTD's mixture fit does not use spatial coordinates for
smoothing; placeholders do not change fitted fractions.

### CIBERSORTx

Obtain access and a token from [CIBERSORTx](https://cibersortx.stanford.edu/), then
set `VHD_CIBERSORTX_USERNAME` and `VHD_CIBERSORTX_TOKEN` in the viewer server's
environment. Do not put credentials into manifests or analysis parameters.
CIBERSORTx requires Docker or Apptainer and uses the official
`cibersortx/fractions` image. For a fixed version, configure a locally approved SIF
or pin the manifest image to the approved digest.

The host exports reference cells using phenotype-labelled columns and invokes the
official single-cell workflow to generate a signature and estimate fractions.
Subsequent layers reuse that signature. Signature generation settings are exposed
in the tab; quantile normalization defaults to false for RNA sequencing. Batch
correction is explicitly selectable (`none`, `S`, or `B`); choose it for the
reference/target technologies in your experiment. Signature and diagnostic files
remain in the run directory. Local Python execution is unavailable for this method.

Commands, captured output, and error tails redact the configured username and
token. Credentials are never stored in AnnData, cache keys, or `params.json`.
Runtime launch arguments necessarily pass credentials to the official executable.
See the [authors' container integration](https://github.com/digitalcytometry/ecotyper)
and the authenticated CIBERSORTx documentation for account/container requirements.

## Service API, cache, and persistence

```python
from hd_viewer.board.operations.layer.deconvolution import DeconvolutionService

service = DeconvolutionService({'runner': 'local', 'work_dir': '/tmp/my-deconvolution'})
result = service.run('nnls', target, reference, 'cell_type', {'batch_size': 256})
# result.fractions, result.genes, result.abundance, result.metadata
```

`run` does not mutate either input. Expression snapshots export only counts,
identifiers, and required reference annotations. Results must exactly match the
exported target observation order and reference cell types. Invalid values,
partial missing rows, duplicate identifiers, or inconsistent abundance/fractions
fail before publication. Zero target rows are omitted from model fitting and
restored as NaN rows. The original `deconvolve(target, reference, label_column,
batch_size=256)` entry point remains available from
`hd_viewer.board.operations.layer.deconvolution`.

In the viewer, runs live in `<dataset base directory>/deconvolution/<method>/<hash>/`,
next to the dataset's `segmentation/` folder. The base directory follows segmentation:
for configured datasets it is the tissue image's parent directory; for uploads it
is the session's uploaded-files directory (`tmp/<session_id>/`). All bins, spots, and segmentation-derived cells
for that dataset use this directory, including layers added later. Reference
artifacts and locks also live under its `deconvolution/` folder. Results remain
stored in the layer's AnnData as well.

Standalone service callers can set constructor `config['work_dir']` or supply
`run(..., work_dir=path)` for a particular dataset. A per-run path does not change
the service's shared queue or another dataset's destination. Standalone/unbound
viewer layers retain the configured fallback directory (`tmp/deconvolution`).
Previously cached runs under the old fallback are not moved automatically. To migrate,
copy the relevant `<method>/<hash>/` folders and `_references/` into the dataset
`deconvolution/` folder, preserving names. Skip `_locks/` and unfinished `.running-*`
folders. Cache loading resolves logs relative to the copied run directory. Verify
that rerunning the same method and parameters reports `Loaded cached result` before
removing the original copy. If the old cache contains multiple datasets, migrate
each dataset's run folders to its own destination.

The content key includes expression, ordered identifiers,
required annotations, selected label column, method version/image, parameters,
and execution identity. Unrelated metadata does not invalidate results. Increment
a plugin's version when changing its implementation. Reference signatures from
cell2location/CIBERSORTx live under `_references/`, independently of spatial data.
Reference artifacts publish only after a complete validated run. File locks avoid
duplicate fitting/publication across processes using the same working directory;
`max_parallel` is a process-local resource budget.

Each successful directory contains `run.log`, result tables, metadata, and any
method artifacts. Failed `.running-*` directories retain inputs and logs for
debugging and are never treated as completed cache entries.

Viewer storage:

- `obs['deconv_<n>:<cell type>']`: fraction columns; dominant labels remain categorical.
- `obs['deconv_<n>_abundance:<position>']`: optional abundance columns, mapped to types in metadata.
- `obsm['deconvolution']`, `uns['deconvolution']`: latest fractions and metadata.
- `obsm['deconvolution_<n>']`: fractions and dominant label for each run.
- `obsm['deconvolution_<n>_abundance']`: optional native abundance for each run.
- `uns['deconvolution_runs'][str(n)]`: method/version, JSON-encoded parameters,
  reference identity, genes, units, log path, and column mappings.

`obsm['deconvolution_abundance']` exists only when the latest run has abundance.
JSON parameter encoding preserves optional null values across AnnData export.

## Creating a plugin

Create a folder containing `manifest.yaml`, `run.py`, `requirements.txt`, and an
optional Dockerfile. The name must match the folder. A minimal manifest:

```yaml
name: example
title: Example method
version: '0.1.0'
image: your-registry/example:0.1.0
description: Explain the estimator and its input requirements.
runners: [docker, apptainer, local]
params:
  - name: batch_size
    type: int
    default: 256
    min: 1
reference_params: []
```

Parameter types are `int`, `float`, `bool`, `choice` (with `choices`), and
`obs_column`. Set `required: true` for required values and `multiple_subjects:
true` on a donor annotation. Optional values default to null. `reference_params`
lists parameters affecting fitted reference artifacts; changing any of them
invalidates the reference cache. Supply `url` and `citation` for the viewer.

The runner mounts read-only inputs at `/data/input` and writable outputs at
`/data/output`; manifests may override `input_mount`, `output_mount`, and
container `command`. Local plugins receive `--input` and `--output` arguments
and the SDK on `PYTHONPATH`. Copy the SDK into container images as in the bundled
Dockerfiles.

Inputs are `target.h5ad`, `reference.h5ad`, `params.json`, and optional read-only
`reference_artifacts/`. Outputs are indexed numeric `fractions.json`, optional
`abundance.json`, and `result.json` containing `genes`, `fraction_units`, and
`abundance_units` when applicable. Use `sdk/vhd_deconv.py` to preserve string
barcodes and labels (including numeric-looking strings) and encode NaN as JSON
null. Plugins must not add the dominant-type column; the viewer derives it.

```python
from vhd_deconv import start
job = start()
target, reference = job.read_inputs()
# Estimate a numeric DataFrame indexed by target.obs_names.
job.report(0.5, 'Fitting')
job.write_result(fractions, genes)
```

Read cached signatures from `job.reference_cache`; write newly fitted reusable
artifacts to `job.reference_output`. They must depend only on reference inputs and
`reference_params`, never on the target. Progress lines use
`##progress <fraction> <message>`. Nonzero exit codes propagate with the log tail.

## Validation

Run from the repository root:

```bash
.venv/bin/python -m unittest tests.test_deconvolution tests.test_segmentation tests.test_deconvolution_plugins tests.test_spacexr_plugin tests.test_param_form
```

The suite runs NNLS as a real isolated subprocess and tests command construction,
CIBERSORTx format conversion, output validation, cache invalidation, reference reuse,
failure recovery, UI controls, history, and AnnData round trips. To run small real
external-method smoke tests, set `VHD_TEST_CELL2LOCATION_PYTHON` and/or
`VHD_TEST_MUSIC_PYTHON` / `VHD_TEST_SPACEXR_PYTHON` to suitable plugin interpreters.
Without those environments, the corresponding tests skip explicitly. MuSiC and
spacexr also need Rscript and their R packages. The spacexr smoke test exercises
all three RCTD modes; contract tests simulate R output and do not verify model fitting.

CIBERSORTx format/command tests use simulated official outputs. Set `VHD_TEST_CIBERSORTX_RUNNER=docker` (or `apptainer`) with server credentials
to run the optional real-container smoke test. A successful credentialed run is required before claiming
the CIBERSORTx integration is verified end to end. Check that the generated
signature and fractions are cached and that a second layer reuses the signature.

Method references: [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html),
[MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html),
[spacexr/RCTD](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html),
[CIBERSORTx](https://cibersortx.stanford.edu/).
