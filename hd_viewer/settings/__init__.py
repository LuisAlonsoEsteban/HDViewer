
from pathlib import Path
import os
from bokeh.palettes import Category20, Category20b, Category20c, Dark2, Muted, Set1


# Paths and services
BASE_DIR = Path(__file__).resolve().parents[1]
DOWNLOADABLE_DATA_PATH = (BASE_DIR / "static" / "data").resolve()
SERVER_HOST = "0.0.0.0"
TILE_SERVER_PORT = 4000
DOWNLOAD_SERVER_PORT = 5002

# Data loading and startup progress
UPLOAD_DIR = (BASE_DIR / ".." / "tmp").resolve()
UPLOAD_BIN_SIZES = (2, 8, 16)
UPLOAD_DEFAULT_BIN_SIZE = 16
PROGRESS_WIDTH = 560
STARTUP_PROGRESS_REFRESH_MS = 50

# Image tiles and grid fitting
TILE_SIZE = 256
TILE_MEM_THRESHOLD = 0.5 # GB
GRID_FIT_MAX_POINTS = 4096
GRID_MAX_CONDITION_NUMBER = 1e10

# Bins and aggregation
DEFAULT_BIN_SIZE = 2.0
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
DEFAULT_CANVAS_PX = 800
VIEWPORT_FETCH_WORKERS = 4
GENE_CACHE_SIZE = 8
assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"

# Clustering, embeddings and expression preprocessing
CLUSTER_COLOR_LABEL = "cluster_color"
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
MIN_GENE_CELLS = 3
NORMALIZATION_TARGET_SUM = 1e4

# Deconvolution
DECONV_RESULT_PREFIX = "deconv_"
DECONV_BATCH_SIZE = 256

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

# Map appearance and panel layout
PLOT_MIN_HEIGHT = 180
PLOT_BACKGROUND_COLOR = "#7e7e7e"
PANEL_SPLITS = {
    "columns": {"sizes": [30, 45, 25], "minSize": 100, "gutterSize": 5},
    "metadata": {"sizes": [50, 50], "minSize": 100, "gutterSize": 5},
    "map": {"sizes": [75, 25], "minSize": [200, 80], "gutterSize": 8},
    "right": {"sizes": [50, 50], "minSize": 100, "gutterSize": 5},
}
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
NOT_SEQUENCED_VAL = {'spots': {'color': '#cccccc'}, 'bins': {'color': '#cccccc'}, 'cells': {'color': '#cccccc'}}
CELL_LAYER_COLORS = ['#6baed6', '#fd8d3c', '#74c476', '#9e9ac8', '#f768a1', '#fdd0a2']
PARAMETER_NOT_SEQUENCED_COLOR = "#f5f5f5"
PARAMETER_COLOR_SPEC = dict(value="#ffffff", not_sequenced="#cccccc")
PARAMETER_ALPHA_SPEC = dict(low=0, high=1, value=1.0, not_sequenced=0.3, minv=0.0, maxv=1.0)
PARAMETER_WIDTH_SPEC = dict(low=0, high=20, value=1.0, not_sequenced=0.5, minv=0.5, maxv=4.0)
CATEGORICAL_PALETTE = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff']
PALETTE_EXCLUDED_GREYS = {'#c7c7c7', '#d9d9d9', '#7f7f7f', '#999999', '#666666', '#bdbdbd', '#969696', '#636363'}
CLUSTER_PALETTE = [c for c in [*Category20[20], *Category20b[20], *Category20c[20]]
                   if c.lower() not in PALETTE_EXCLUDED_GREYS]
GROUP_PALETTE = [c for c in [*Set1[9], *Dark2[8], *Muted[9]]
                 if c.lower() not in PALETTE_EXCLUDED_GREYS]

# Text autocomplete
AUTOCOMPLETE_MIN_CHARACTERS = 2
METADATA_INPUT_MIN_CHARACTERS = 1
AUTOCOMPLETE_CASE_SENSITIVE = False
AUTOCOMPLETE_SEARCH_STRATEGY = "starts_with"

# Loading status
MAX_LOG_LINES = 6
LOADING_FONT_SIZE = 12
LOADING_REFRESH_MS = 500
