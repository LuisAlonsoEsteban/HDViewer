"""Discovery of segmentation plugins (segmentation/plugins/<name>/manifest.yaml)."""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml

SEGMENTATION_ROOT = Path(__file__).resolve().parents[2] / "segmentation"
PLUGINS_DIR = SEGMENTATION_ROOT / "plugins"
SDK_DIR = SEGMENTATION_ROOT / "sdk"

PARAM_TYPES = ("float", "int", "bool", "choice")


@dataclass
class Param:
    name: str
    type: str
    default: object = None
    label: str = None
    help: str = None
    min: float = None
    max: float = None
    step: float = None
    choices: list = None

    def __post_init__(self):
        if self.type not in PARAM_TYPES:
            raise ValueError(f"parameter '{self.name}': type must be one of {PARAM_TYPES}")
        if self.type == "choice" and not self.choices:
            raise ValueError(f"parameter '{self.name}': a choice needs 'choices'")
        self.label = self.label or self.name

    def coerce(self, value):
        if value is None or value == "":
            return None
        if self.type == "float":
            value = float(value)
        elif self.type == "int":
            value = int(value)
        elif self.type == "bool":
            value = bool(value)
        elif value not in self.choices:
            raise ValueError(f"'{value}' is not a valid {self.label}: choose from {self.choices}")
        if self.type in ("float", "int"):
            if self.min is not None and value < self.min:
                raise ValueError(f"{self.label} must be ≥ {self.min}")
            if self.max is not None and value > self.max:
                raise ValueError(f"{self.label} must be ≤ {self.max}")
        return value


@dataclass
class Plugin:
    name: str
    title: str
    path: Path
    image: str
    version: str = "0"
    description: str = ""
    url: str = None
    citation: str = None
    params: list = field(default_factory=list)

    def resolve_params(self, overrides=None):
        """Defaults from the manifest, overridden and validated by `overrides`."""
        overrides = overrides or {}
        unknown = set(overrides) - {p.name for p in self.params}
        if unknown:
            raise ValueError(f"{self.title}: unknown parameters {sorted(unknown)}")
        return {p.name: p.coerce(overrides.get(p.name, p.default)) for p in self.params}


def load_plugin(path):
    path = Path(path)
    with open(path / "manifest.yaml") as f:
        manifest = yaml.safe_load(f) or {}

    for key in ("name", "title", "image"):
        if not manifest.get(key):
            raise ValueError(f"{path / 'manifest.yaml'}: '{key}' is required")
    if manifest["name"] != path.name:
        raise ValueError(f"{path / 'manifest.yaml'}: name '{manifest['name']}' must match the folder name")
    if not (path / "run.py").exists():
        raise ValueError(f"{path}: run.py is missing")

    params = [Param(**spec) for spec in manifest.get("params") or []]
    return Plugin(
        name=manifest["name"], title=manifest["title"], path=path,
        image=manifest["image"], version=str(manifest.get("version", "0")),
        description=manifest.get("description", ""), url=manifest.get("url"),
        citation=manifest.get("citation"), params=params,
    )


def discover(plugins_dir=PLUGINS_DIR):
    """All valid plugins by name; invalid ones are reported and skipped."""
    plugins = {}
    for path in sorted(Path(plugins_dir).iterdir()):
        if not (path / "manifest.yaml").exists():
            continue
        try:
            plugin = load_plugin(path)
        except Exception as exc:
            print(f"[segmentation] skipping plugin {path.name}: {exc}", file=sys.stderr)
            continue
        plugins[plugin.name] = plugin
    return plugins
