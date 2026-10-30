
from pathlib import Path
from .....plugins.manifest import PARAM_TYPES, Param, Plugin, load_plugin, discover as _discover

SEGMENTATION_ROOT = Path(__file__).resolve().parent / "api"
PLUGINS_DIR = SEGMENTATION_ROOT / "plugins"
SDK_DIR = SEGMENTATION_ROOT / "sdk"


def discover(plugins_dir=PLUGINS_DIR):
    return _discover(plugins_dir, kind="segmentation")
