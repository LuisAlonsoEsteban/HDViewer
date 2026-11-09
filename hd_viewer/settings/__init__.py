
from pathlib import Path
import os


# Paths and services
BASE_DIR = Path(__file__).resolve().parents[1]
# Root of the bundled example datasets (see hd_viewer/startup/datasets.py).
DATA_DIR = Path(os.environ.get("VHD_DATA_DIR", BASE_DIR.parent / "data")).resolve()

# Data loading and startup progress
UPLOAD_DIR = (BASE_DIR / ".." / "tmp").resolve()
UPLOAD_BIN_SIZES = (2, 8, 16)
UPLOAD_DEFAULT_BIN_SIZE = 16
# Largest single file accepted by the startup uploads (bytes; default 20 GiB).
# Files are streamed to disk in chunks, so this bounds disk use, not memory.
UPLOAD_MAX_BYTES = int(os.environ.get("VHD_UPLOAD_MAX_BYTES", 20 * 1024 ** 3))
# Size of each upload chunk (one HTTP PUT of raw bytes); keep it below the reverse
# proxy's request body limit (deploy/nginx-hdviewer.conf allows 16 MiB).
UPLOAD_CHUNK_BYTES = 4 * 1024 ** 2

# Image tiles and grid fitting
TILE_SIZE = 256
TILE_MEM_THRESHOLD = 0.5 # GB
# Tile encoding: "webp" (lossy, ~15x smaller than PNG) or "png" (lossless).
# Changing either value regenerates the tile caches (they are part of the manifest).
TILE_FORMAT = os.environ.get("VHD_TILE_FORMAT", "webp").lower()
TILE_QUALITY = int(os.environ.get("VHD_TILE_QUALITY", 80))  # WebP quality 0-100; ignored for PNG
GRID_FIT_MAX_POINTS = 4096
GRID_MAX_CONDITION_NUMBER = 1e10

# Bins and aggregation
AGGREGATED_PREFIX = "agg:"
AGGREGATE_OBS_PREFIX = "summary:"
AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
MAX_CACHE_ENTRIES = 24

# Spot buffer detail
SPOT_BUFFER_QUAD_SEGS = 16

# Cell levels of detail
LOD_PREFIX = "lod:"
LOD_BASE_TOLERANCE = 1.0
LOD_MAX_LEVEL = 8
LOD_TOLERANCE_PX = 1.0
LOD_HYSTERESIS = 0.15

# Viewport loading and performance
STATIC_RANGE = 0.6
OFFSET = 0.7
HYSTERESIS_RATIO = 0.15
TARGET_BINS_PER_AXIS = 180
GENE_CACHE_SIZE = 8
# Far rows are dropped only once the loaded rows exceed this multiple of the
# last viewport fetch (each prune resends every retained row)
PRUNE_ROW_BUDGET = 3
assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"

# Clustering, embeddings and expression preprocessing
CLUSTER_RESOLUTION = 0.9
CLUSTER_ALGORITHM = "leiden"
CLUSTER_N_CLUSTERS = 8
CLUSTER_DBSCAN_EPS = 0.5
CLUSTER_MIN_SAMPLES = 5
CLUSTER_MIN_SIZE = 5
CLUSTER_COVARIANCE_TYPE = "full"
CLUSTER_LINKAGE = "ward"
CLUSTER_RANDOM_STATE = 0
CLUSTER_SLIDER_DEFAULT = 1.0
CLUSTER_SLIDER_MIN = 0.1
CLUSTER_SLIDER_MAX = 2.0
CLUSTER_SLIDER_STEP = 0.1
N_PCS = 30
N_NEIGHBORS = 15
MARKER_GENE_COUNT = 5
METADATA_MARKER_GENE_COUNT = 10
QC_SCATTER_MAX_POINTS = 50000        # React UI: QC scatters show a random subset beyond this
MIN_GENE_CELLS = 3
NORMALIZATION_TARGET_SUM = 1e4

# Deconvolution
DECONV_RESULT_PREFIX = "deconv_"

# Groups and differential expression
GROUPS_COL = "custom_groups"
GROUP_MASK_PREFIX = "grp_"
DE_CLUSTER_PREFIX = "Cluster "
DE_METHOD = "wilcoxon"
DE_MIN_OBSERVATIONS = 3
DE_PVALUE_CUTOFF = 0.05
DE_TABLE_ROWS = 200
DE_PLOT_TOP_GENES = 5
DE_KDE_MAX_SAMPLES = 5000
DE_RANDOM_SEED = 0
DE_VIOLIN_GRID_POINTS = 128
DE_VIOLIN_HALF_WIDTH = 0.28

# Sample comparison
COMPARISON_OBS_PREFIX = "cmp_"
COMPARISON_PREDICATE = "intersects"
COMPARISON_MIN_COVERAGE = 0.8
COMPARISON_COVERAGE_STEP = 0.05
COMPARISON_N_HVG = 2000
COMPARISON_TOP_SPATIAL_GENES = 100
COMPARISON_CHUNK_SIZE = 256
COMPARISON_K_NEIGHBORS = 6
COMPARISON_HULL_RATIO = 0.1

# ROI interaction and background processing
ROI_DEFAULT_RADIUS = 50.0
ROI_FOCUS_PADDING = 0.05
ROI_BUFFER_QUAD_SEGS = 64
ROI_SEED_BATCH_SIZE = 256
ROI_PREVIEW_DEBOUNCE_MS = 200
ROI_BARCODE_MAX_COMPLETIONS = 50
ROI_OVERLAY_MAX_VERTICES = 20000
ROI_OVERLAY_MAX_SIZE = 768

# Segmentation; config.yaml and explicit parameters override these defaults
SEGMENTATION_ROOT = (BASE_DIR / "startup" / "data" / "image" / "segmentation" / "api").resolve()
SEGMENTATION_CONFIG_PATH = SEGMENTATION_ROOT / "config.yaml"
PLUGINS_DIR = SEGMENTATION_ROOT / "plugins"
SDK_DIR = SEGMENTATION_ROOT / "sdk"
SEGMENTATION_DEFAULT_RUNNER = "auto"
DEFAULT_THREADS = 8
SEGMENTATION_MAX_PARALLEL = 2
SEGMENTATION_GPU = False
SEGMENTATION_DOCKER_BINARY = "docker"
SEGMENTATION_APPTAINER_BINARY = "apptainer"
SEGMENTATION_IMAGES_DIR = "images"
# Per run: wall-clock limit (minutes), and container RAM / process caps (Docker)
SEGMENTATION_TIMEOUT_MINUTES = 120
SEGMENTATION_MEMORY = "32g"
SEGMENTATION_PIDS_LIMIT = 4096

# Map appearance and panel layout
PLOT_MIN_HEIGHT = 180
PLOT_BACKGROUND_COLOR = "#0f1d17"
RULER_COLOR = "#ffd166"

# Layer colours and drawing parameter defaults
INITIAL_ATTRS = {
    "color": "#6baed6",
    "alpha": 0.5,
}
NOT_SEQUENCED_ATTRS = {
    "color": "#b8b8b8",
    "alpha": 0.5,
}
CELL_LAYER_COLORS = ['#6baed6', '#fd8d3c', '#74c476', '#9e9ac8', '#f768a1', '#fdd0a2']
PARAMETER_NOT_SEQUENCED_COLOR = "#f5f5f5"
PARAMETER_COLOR_SPEC = dict(value="#ffffff", not_sequenced="#cccccc")
PARAMETER_ALPHA_SPEC = dict(low=0, high=1, value=1.0, not_sequenced=0.3, minv=0.0, maxv=1.0)
PARAMETER_WIDTH_SPEC = dict(low=0, high=20, value=1.0, not_sequenced=0.5, minv=0.5, maxv=4.0)
CATEGORICAL_PALETTE = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff']
PALETTE_EXCLUDED_GREYS = {'#c7c7c7', '#d9d9d9', '#7f7f7f', '#999999', '#666666', '#bdbdbd', '#969696', '#636363'}
# Bokeh's Category20/20b/20c, Set1, Dark2 and Muted palettes, copied so settings
# (and every backend module) import without Bokeh
_CATEGORY20 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
               '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
_CATEGORY20B = ['#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39',
                '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b', '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6']
_CATEGORY20C = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2', '#31a354', '#74c476',
                '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc', '#dadaeb', '#636363', '#969696', '#bdbdbd', '#d9d9d9']
_SET1 = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
_DARK2 = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']
_MUTED = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499']
CLUSTER_PALETTE = [c for c in [*_CATEGORY20, *_CATEGORY20B, *_CATEGORY20C]
                   if c.lower() not in PALETTE_EXCLUDED_GREYS]
GROUP_PALETTE = [c for c in [*_SET1, *_DARK2, *_MUTED]
                 if c.lower() not in PALETTE_EXCLUDED_GREYS]


# Deconvolution plugins have independent execution and cache configuration.
DECONVOLUTION_ROOT = BASE_DIR / 'analysis' / 'deconvolution' / 'api'
DECONVOLUTION_CONFIG_PATH = DECONVOLUTION_ROOT / 'config.yaml'
# Fallback for standalone layers; loaded datasets supply <data path>/deconvolution.
DECONVOLUTION_WORK_DIR = UPLOAD_DIR / 'deconvolution'
DECONVOLUTION_DEFAULT_RUNNER = 'auto'
DECONVOLUTION_GPU = False
# RCTD on a full Visium HD slide takes ~8 h
DECONVOLUTION_TIMEOUT_MINUTES = 720
DECONVOLUTION_MEMORY = "64g"
DECONVOLUTION_PIDS_LIMIT = 4096
DECONVOLUTION_MAX_PARALLEL = 1
