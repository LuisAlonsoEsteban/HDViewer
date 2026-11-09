"""Discover independently packaged deconvolution methods."""
from pathlib import Path
from ...plugins.manifest import Param, Plugin, load_plugin, discover as _discover

ROOT = Path(__file__).resolve().parent / 'api'
PLUGINS_DIR = ROOT / 'plugins'
SDK_DIR = ROOT / 'sdk'


def discover(plugins_dir=PLUGINS_DIR):
    return _discover(plugins_dir, kind='deconvolution')
