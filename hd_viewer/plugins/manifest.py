"""Manifest discovery shared by analysis plugins."""


import sys
import math
from dataclasses import dataclass, field
from pathlib import Path

import yaml

PARAM_TYPES = ('float', 'int', 'bool', 'choice', 'obs_column')


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
    required: bool = False
    multiple_subjects: bool = False

    def __post_init__(self):
        if self.type not in PARAM_TYPES:
            raise ValueError(f"parameter '{self.name}': type must be one of {PARAM_TYPES}")
        if self.type == "choice" and not self.choices:
            raise ValueError(f"parameter '{self.name}': a choice needs 'choices'")
        self.label = self.label or self.name

    def coerce(self, value):
        if value is None or value == "":
            if self.default is not None:
                value = self.default
            elif self.required:
                raise ValueError(f"{self.label} is required")
            else:
                return None
        if self.type == "float":
            value = float(value)
        elif self.type == "int":
            number = float(value)
            if not math.isfinite(number) or not number.is_integer():
                raise ValueError(f"{self.label} must be an integer")
            value = int(number)
        elif self.type == "bool":
            if isinstance(value, str):
                if value.lower() not in ("true", "false"):
                    raise ValueError(f"{self.label} must be true or false")
                value = value.lower() == "true"
            elif not isinstance(value, bool):
                raise ValueError(f"{self.label} must be true or false")
        elif self.type == "obs_column":
            value = str(value)
        elif value not in self.choices:
            raise ValueError(f"'{value}' is not a valid {self.label}: choose from {self.choices}")
        if self.type in ("float", "int"):
            if not math.isfinite(value):
                raise ValueError(f"{self.label} must be finite")
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

    runners: list = field(default_factory=lambda: ["docker", "apptainer", "local"])
    command: list = field(default_factory=list)
    input_mount: str = "/data/input"
    output_mount: str = "/data/output"
    reference_params: list = field(default_factory=list)
    # flag -> server environment variable whose value the runner appends as
    # `flag value` without putting it in any command line it launches
    secret_args: dict = field(default_factory=dict)

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
    if len({p.name for p in params}) != len(params):
        raise ValueError("Duplicate parameter names")
    if set(manifest.get("reference_params", [])) - {p.name for p in params}:
        raise ValueError("reference_params must name declared parameters")
    runners = manifest.get("runners", ["docker", "apptainer", "local"])
    if not runners or set(runners) - {"docker", "apptainer", "local"}:
        raise ValueError("Invalid supported runners")
    
    return Plugin(
        name=manifest["name"], title=manifest["title"], path=path,
        image=manifest["image"], version=str(manifest.get("version", "0")),
        description=manifest.get("description", ""), url=manifest.get("url"),
        citation=manifest.get("citation"), params=params,
        runners=runners, command=manifest.get("command", []),
        input_mount=manifest.get("input_mount", "/data/input"),
        output_mount=manifest.get("output_mount", "/data/output"),
        reference_params=manifest.get("reference_params", []),
    )


def discover(plugins_dir, *, kind="plugin"):
    """All valid plugins by name; invalid ones are reported and skipped."""
    plugins = {}
    if not Path(plugins_dir).is_dir():
        return plugins
    for path in sorted(Path(plugins_dir).iterdir()):
        if not (path / "manifest.yaml").exists():
            continue
        try:
            plugin = load_plugin(path)
        except Exception as exc:
            print(f"[{kind}] skipping plugin {path.name}: {exc}", file=sys.stderr)
            continue
        plugins[plugin.name] = plugin
    return plugins
