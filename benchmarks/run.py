"""Run the cross-tool benchmark (run with the project's .venv: needs psutil, pyyaml, h5py).

    .venv/bin/python benchmarks/run.py                          # everything in config.yaml
    .venv/bin/python benchmarks/run.py --datasets syn_100k --repeats 1
    .venv/bin/python benchmarks/run.py --tools hdviewer,spatialdata --tasks 'bins.*' --datasets mouse_brain
    .venv/bin/python benchmarks/run.py --list                   # show the planned runs only

Each (tool, dataset, task, bin size, repeat) runs in a fresh process with the
tool's interpreter (harness/measure.py), and appends one JSON line to
results/<name>.jsonl. Summaries and figures: benchmarks/report.py.
"""
import argparse
import fnmatch
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

BENCH = Path(__file__).resolve().parent
REPO = BENCH.parent
sys.path.insert(0, str(REPO))

from benchmarks.harness import datasets as ds  # noqa: E402
from benchmarks.harness import synthetic  # noqa: E402
from benchmarks.harness.measure import run_task  # noqa: E402


def _path(p):
    p = Path(p)
    return p if p.is_absolute() else REPO / p


def synthetic_entries(cfg):
    syn = cfg.get("synthetic") or {}
    out = {}
    for size in syn.get("sizes", []):
        root = Path(syn["out"]) / str(size)
        out[f"syn_{size}"] = dict(root=str(root), image=str(root / "synthetic_tissue_image.tif"),
                                  cells=str(root / "cells.geojson"), dataset_id="synthetic",
                                  genes=syn["marker_genes"], bin_sizes=syn["bin_sizes"],
                                  n_bins_2um=synthetic.parse_size(size), synthetic=dict(size=size, **syn))
    return out


def ensure_synthetic(entry):
    root = _path(entry["root"])
    if (root / "meta.json").exists():
        return
    syn = entry["synthetic"]
    print(f"[synthetic] generating {syn['size']} bins in {root}", flush=True)
    synthetic.generate(synthetic.parse_size(syn["size"]), root, n_genes=syn["genes"], seed=syn["seed"])


def adapter_for(cfg, tool, task):
    """The adapter listing `task` explicitly, else the tool's general adapter."""
    adapters = cfg["tools"][tool].get("adapters", [])
    for adapter in adapters:
        if task in adapter.get("tasks", ()):
            return adapter
    for adapter in adapters:
        if "tasks" not in adapter:
            return adapter
    return None


def all_tasks(cfg):
    return [t for group in cfg["tasks"].values() for t in group]


def bin_sizes_for(cfg, task, spec):
    if not task.startswith(("bins.", "e2e.")):
        return [None]
    wanted = cfg.get("task_bin_sizes", {}).get(task, spec["bin_sizes"])
    return [b for b in wanted if b in spec["bin_sizes"]]


class Xvfb:
    """One X server (OpenGL through Mesa) for the napari runs."""

    def __init__(self, cfg):
        self.cfg = cfg.get("xvfb") or {}
        self.proc = None

    def env(self):
        if self.proc is None:
            for display in range(90, 140):
                if not Path(f"/tmp/.X11-unix/X{display}").exists() and not Path(f"/tmp/.X{display}-lock").exists():
                    break
            prefix = Path(self.cfg["prefix"])
            env = dict(os.environ, LIBGL_DRIVERS_PATH=self.cfg.get("dri_path", "/usr/lib64/dri"))
            self.proc = subprocess.Popen([str(prefix / "bin" / "Xvfb"), f":{display}", "-screen", "0",
                                          "1920x1200x24", "+extension", "GLX", "-nolisten", "tcp"],
                                         env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            time.sleep(2)
            self.display = f":{display}"
        prefix = Path(self.cfg["prefix"])
        return dict(DISPLAY=self.display, QT_API="pyside6", QT_QPA_PLATFORM="xcb",
                    LD_LIBRARY_PATH=str(prefix / "qtlibs"))

    def close(self):
        if self.proc is not None:
            self.proc.terminate()


def machine():
    import psutil
    return dict(host=socket.gethostname(), cpus=os.cpu_count(), memory_gb=round(psutil.virtual_memory().total / 2**30),
                platform=platform.platform(), git=subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                                                                 capture_output=True, text=True).stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(BENCH / "config.yaml"))
    parser.add_argument("--tools", help="comma list (default: all automated tools)")
    parser.add_argument("--datasets", help="comma list: mouse_brain, syn_100k, ... (default: all)")
    parser.add_argument("--tasks", help="comma list of task names or globs, e.g. 'bins.*,image.render'")
    parser.add_argument("--bin-sizes", help="comma list restricting bins/e2e tasks, e.g. 8,16")
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--timeout", type=int)
    parser.add_argument("--name", help="results file name (default: date-time)")
    parser.add_argument("--list", action="store_true", help="print the planned runs and exit")
    parser.add_argument("--resume", action="store_true",
                        help="skip runs already recorded in results/<name>.jsonl (needs --name)")
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config))
    repeats = args.repeats or cfg["repeats"]
    timeout = args.timeout or cfg["timeout_s"]
    entries = dict(cfg.get("datasets") or {}, **synthetic_entries(cfg))
    names = args.datasets.split(",") if args.datasets else list(entries)
    tools = args.tools.split(",") if args.tools else [t for t, c in cfg["tools"].items() if not c.get("manual")]
    patterns = args.tasks.split(",") if args.tasks else ["*"]
    tasks = [t for t in all_tasks(cfg) if any(fnmatch.fnmatch(t, p) for p in patterns)]
    only_bins = [int(b) for b in args.bin_sizes.split(",")] if args.bin_sizes else None

    work = _path(cfg["work_dir"])
    cache = _path(cfg["cache_dir"])
    results = _path(cfg["results_dir"])
    results.mkdir(parents=True, exist_ok=True)
    out = results / f"{args.name or datetime.now().strftime('%Y%m%d-%H%M%S')}.jsonl"

    plan = []
    specs = {}
    for name in names:
        if name not in entries:
            sys.exit(f"Unknown dataset {name}; choose from {', '.join(entries)}")
        entry = entries[name]
        if not args.list and "synthetic" in entry:
            ensure_synthetic(entry)
        if args.list:
            specs[name] = dict(bin_sizes=entry["bin_sizes"])
        else:
            specs[name] = ds.prepare(name, entry, work, cfg["analysis"])
        for task in tasks:
            if task.startswith("cells.") and not entry.get("cells"):
                continue
            for tool in tools:
                adapter = adapter_for(cfg, tool, task)
                if adapter is None:
                    continue
                for b in bin_sizes_for(cfg, task, specs[name]):
                    if only_bins and b is not None and b not in only_bins:
                        continue
                    plan.append((name, task, b, tool, adapter))
    print(f"{len(plan)} runs x {repeats} repeats -> {out}", flush=True)
    if args.list:
        for name, task, b, tool, _ in plan:
            print(f"  {name:14s} {task:16s} {'' if b is None else f'{b} µm':6s} {tool}")
        return

    threads = str(cfg["threads"])
    base_env = {k: str(v) for k, v in (cfg.get("env") or {}).items()}
    base_env.update(OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads, OPENBLAS_NUM_THREADS=threads,
                    NUMBA_NUM_THREADS=threads, MPLBACKEND="Agg", PYTHONWARNINGS="ignore")
    host = machine()
    xvfb = Xvfb(cfg)
    try:
        done = set()
        if args.resume and out.exists():
            for line in open(out):
                try:
                    r = json.loads(line)
                except ValueError:
                    continue        # a line cut by a killed job
                done.add((r["tool"], r["dataset"], r["task"], r["bin_size"], r["repeat"]))
                if r["status"] == "n.a.":
                    done.update((r["tool"], r["dataset"], r["task"], r["bin_size"], k) for k in range(repeats))
            print(f"resuming: {len(done)} runs already recorded", flush=True)
        with open(out, "a") as sink:
            for i, (name, task, b, tool, adapter) in enumerate(plan, 1):
                for repeat in range(repeats):
                    if (tool, name, task, b, repeat) in done:
                        continue
                    run_dir = work / "runs" / f"{tool}-{name}-{task}-{b}-{repeat}"
                    shutil.rmtree(run_dir, ignore_errors=True)
                    run_dir.mkdir(parents=True)
                    spec = dict(specs[name], bin_size=b or specs[name]["bin_sizes"][-1], workdir=str(run_dir),
                                cache_dir=str(cache), hdviewer_python=str(_path(cfg["python"]["hdviewer"])),
                                snapshots=cfg.get("snapshots", False))
                    env = dict(base_env, **(xvfb.env() if adapter.get("xvfb") else {}))
                    cmd = [str(_path(cfg["python"][adapter["python"]])), str(BENCH / "adapters" / adapter["script"]),
                           "--task", task, "--spec", json.dumps(spec)]
                    result = run_task(cmd, env=env, cwd=str(REPO), timeout=timeout,
                                      cap_bytes=cfg["memory_cap_gb"] * 2**30, hz=cfg["sample_hz"],
                                      log_path=str(run_dir / "run.log"))
                    record = dict(tool=tool, tool_label=cfg["tools"][tool]["label"], dataset=name,
                                  n_bins_2um=specs[name].get("n_bins_2um"), task=task, modality=task.split(".")[0],
                                  bin_size=b, repeat=repeat, time=datetime.now().isoformat(timespec="seconds"),
                                  machine=host, **{k: v for k, v in result.items() if k != "task"})
                    sink.write(json.dumps(record, default=str) + "\n")
                    sink.flush()
                    total = sum(p["seconds"] for p in result.get("phases", []) if not p["name"].startswith("setup:"))
                    print(f"[{i}/{len(plan)} r{repeat}] {tool:11s} {name:12s} {task:15s} "
                          f"{'' if b is None else f'{b}µm':5s} {result['status']:7s} {total:9.2f}s "
                          f"peak {(result.get('peak_tree_bytes') or 0) / 2**30:6.2f} GiB "
                          f"{(result.get('reason') or '')[:80]}", flush=True)
                    if result["status"] == "n.a.":
                        break               # no repeats of a task the tool does not have
                    keep = run_dir / "run.log"
                    for item in run_dir.iterdir():   # free disk: keep logs and pictures only
                        if item != keep and item.suffix not in (".png", ".log"):
                            shutil.rmtree(item, ignore_errors=True) if item.is_dir() else item.unlink()
    finally:
        xvfb.close()


if __name__ == "__main__":
    main()
