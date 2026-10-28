
# Clustering
CLUSTER_PREFIX = "leiden_"
CLUSTER_COL = "leiden"

# Deconvolution
DECONV_RESULT_PREFIX = "deconv_"

# Tiles
TILE_MEM_THRESHOLD = 0.5 # GB

# Bins
AGGREGATED_PREFIX = "agg:"
AGGREGATION_FACTORS = (1, 2, 4, 8, 12)
MAX_CACHE_ENTRIES = 24

# Cells
LOD_PREFIX = "lod:"
LOD_BASE_TOLERANCE = 1.0
LOD_MAX_LEVEL = 8
LOD_TOLERANCE_PX = 1.0
LOD_HYSTERESIS = 0.15

# Lazy loader
STATIC_RANGE = 0.6
OFFSET = 0.7
assert STATIC_RANGE <= OFFSET, "Offset must be larger than the static range for efficient loading"
HYSTERESIS_RATIO = 0.15

# Segmentation
SEGMENTATION_ROOT = Path(__file__).resolve().parents[1] / "segmentation"
PLUGINS_DIR = SEGMENTATION_ROOT / "plugins"
SDK_DIR = SEGMENTATION_ROOT / "sdk"
PARAM_TYPES = ("float", "int", "bool", "choice")

# Load data
UPLOAD_DIR = Path("./tmp")
PROGRESS_WIDTH = 560  # px; fixed so the bars do not resize as their text changes

# Attributes
INITIAL_ATTRS = {
    "color": "#6baed6",
    "alpha": 0.5,
}
NOT_SEQUENCED_ATTRS = {
    "color": "#b8b8b8",
    "alpha": 0.5,
}
CLUSTER_COLOR_LABEL = "cluster_color"

# Downloads
STATIC_PATH = Path("/bokeh_server/static/data")
BASE_DIR = Path(__file__).resolve().parent  
DOWNLOADABLE_DATA_PATH = (BASE_DIR.parent / "static" / "data").resolve()
DOWNLOADABLE_DATA_PATH.mkdir(parents=True, exist_ok=True)

# Groups
GROUPS_COL = "custom_groups"

# Loading
MAX_LOG_LINES = 6