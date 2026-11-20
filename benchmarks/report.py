"""Summaries and figures from benchmark results (run with the project's .venv).

    .venv/bin/python benchmarks/report.py                         # every results/*.jsonl
    .venv/bin/python benchmarks/report.py results/run1.jsonl --loupe adapters/loupe/results.csv

Writes to results/report/:
  phases.csv    one row per tool x dataset x task x bin size x phase: median, IQR, first run, n
  tasks.csv     the same per task (time = sum of non-setup phases; memory = peak of the process tree)
  outcomes.csv  status per tool x dataset x task (ok / n.a. / fail / timeout / oom) with reasons
  summary.md    the key metrics as Markdown tables (paper supplement)
  fig_<dataset>_time.(pdf|png), fig_<dataset>_memory.(pdf|png)   key metrics per tool
  fig_scaling_time.(pdf|png), fig_scaling_memory.(pdf|png)     synthetic sweep (log-log)
  fig_outcomes.(pdf|png)                                         what each tool can do

Setup phases (building caches the task starts from) are excluded from task
figures. The first repeat is reported separately ("first", OS caches cold for
that tool) and the rest as median and IQR.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

BENCH = Path(__file__).resolve().parent

# Categorical slots in fixed order (validated: benchmarks/README.md, "Figures")
TOOLS = ["hdviewer", "spatialdata", "napari", "vitessce", "loupe"]
LABELS = {"hdviewer": "HDViewer (ours)", "spatialdata": "SpatialData", "napari": "napari",
          "vitessce": "Vitessce", "loupe": "Loupe Browser"}
COLOURS = {"hdviewer": "#2a78d6", "spatialdata": "#eb6834", "napari": "#1baf7a", "vitessce": "#eda100",
           "loupe": "#e87ba4"}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

# (title, task, phase or None for the task total, bin size or None, divide by)
KEY_METRICS = [
    ("Image: read", "image.read", None, None, 1),
    ("Image: build pyramid", "image.pyramid", None, None, 1),
    ("Image: overview render", "image.render", "view:overview", None, 1),
    ("Image: pan step", "image.render", "pan", None, 10),
    ("Bins 8 µm: load", "bins.load", None, 8, 1),
    ("Bins 8 µm: viewport query", "bins.query", "view:overview", 8, 1),
    ("Bins 8 µm: overview render", "bins.render", "view:overview", 8, 1),
    ("Bins 8 µm: pan step", "bins.render", "pan", 8, 10),
    ("Bins 8 µm: gene switch", "bins.render", "gene_switch", 8, 1),
    ("Bins 2 µm: load", "bins.load", None, 2, 1),
    ("Bins 2 µm: overview render", "bins.render", "view:overview", 2, 1),
    ("Bins 2 µm: pan step", "bins.render", "pan", 2, 10),
    ("Bins 2 → 8/16 µm: aggregate", "bins.aggregate", None, 2, 1),
    ("Bins 8 µm: clustering + UMAP", "bins.cluster", None, 8, 1),
    ("Bins 8 µm: DE (Wilcoxon)", "bins.de", None, 8, 1),
    ("Cells: load outlines", "cells.load", None, None, 1),
    ("Cells: pool 2 µm bins", "cells.pool", None, None, 1),
    ("Cells: overview render", "cells.render", "view:overview", None, 1),
    ("Cells: pan step", "cells.render", "pan", None, 10),
    ("Cells: clustering + UMAP", "cells.cluster", None, None, 1),
    ("Cells: DE (Wilcoxon)", "cells.de", None, None, 1),
    ("Launch → first view (8 µm)", "e2e.first_view", "load_to_first_view", 8, 1),
]
SCALING = [
    ("Bins 2 µm: load", "bins.load", None, 2, 1),
    ("Bins 2 µm: overview render", "bins.render", "view:overview", 2, 1),
    ("Bins 2 µm: pan step", "bins.render", "pan", 2, 10),
    ("Bins 2 µm: viewport query", "bins.query", "view:overview", 2, 1),
    ("Cells: pool 2 µm bins", "cells.pool", None, None, 1),
    ("Cells: overview render", "cells.render", "view:overview", None, 1),
]


def load(paths, loupe_csv=None):
    records = []
    for path in paths:
        with open(path) as f:
            records += [json.loads(line) for line in f if line.strip()]
    if loupe_csv:
        records += read_loupe(loupe_csv)
    return records


def read_loupe(path):
    """Manual Loupe measurements (adapters/loupe/results_template.csv) as records."""
    grouped = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            if not row.get("task") or row["task"].startswith("#"):
                continue
            key = (row["dataset"], row["task"], row.get("bin_size") or "", row.get("repeat") or "0",
                   row.get("status") or "ok")
            grouped[key].append(row)
    out = []
    for (dataset, task, bin_size, repeat, status), rows in grouped.items():
        phases = [dict(name=r["phase"] or "total", seconds=float(r["seconds"]) if r["seconds"] else None,
                       peak_tree_bytes=float(r["peak_rss_mb"]) * 2**20 if r.get("peak_rss_mb") else None)
                  for r in rows if status == "ok"]
        out.append(dict(tool="loupe", dataset=dataset, task=task, bin_size=int(bin_size) if bin_size else None,
                        repeat=int(repeat), status=status, reason="; ".join(r.get("notes", "") for r in rows),
                        phases=phases, memory_metric="rss (manual)"))
    return out


def task_value(record, phase=None, divide=1):
    """(seconds, peak bytes) of one record: a phase, or the task (non-setup phases)."""
    phases = [p for p in record.get("phases", []) if not p["name"].startswith("setup:")]
    if phase is not None:
        phases = [p for p in phases if p["name"] == phase]
        if not phases:
            return None, None
    if not phases or any(p.get("seconds") is None for p in phases):
        return None, None
    seconds = sum(p["seconds"] for p in phases) / divide
    peaks = [p.get("peak_tree_bytes") for p in phases if p.get("peak_tree_bytes")]
    return seconds, (max(peaks) if peaks else None)


def stats(values):
    """dict(first, median, q1, q3, n) of the values in repeat order."""
    values = [v for v in values if v is not None]
    if not values:
        return None
    rest = values[1:] or values
    q1, med, q3 = np.percentile(rest, [25, 50, 75])
    return dict(first=values[0], median=float(med), q1=float(q1), q3=float(q3), n=len(values))


def group(records):
    out = defaultdict(list)
    for r in sorted(records, key=lambda r: r.get("repeat", 0)):
        out[(r["tool"], r["dataset"], r["task"], r.get("bin_size"))].append(r)
    return out


def outcome(runs):
    statuses = [r["status"] for r in runs]
    for s in ("ok", "timeout", "oom", "fail", "n.a."):
        if s in statuses:
            if s == "ok" and len(set(statuses)) > 1:
                return "partial", "; ".join(sorted({r.get("reason", "") for r in runs if r["status"] != "ok"}))
            return s, (runs[statuses.index(s)].get("reason") or "")
    return statuses[0], ""


def metric(groups, tool, dataset, task, phase, bin_size, divide):
    runs = groups.get((tool, dataset, task, bin_size), [])
    if not runs:
        return None, None, "missing"
    status, _ = outcome(runs)
    ok = [r for r in runs if r["status"] == "ok"]
    values = [task_value(r, phase, divide) for r in ok]
    return stats([v[0] for v in values]), stats([v[1] for v in values]), status


# -- tables ----------------------------------------------------------------------------

def write_tables(groups, out):
    with open(out / "phases.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tool", "dataset", "task", "bin_size", "phase", "seconds_median", "seconds_q1", "seconds_q3",
                    "seconds_first", "peak_gib_median", "n"])
        for (tool, dataset, task, b), runs in sorted(groups.items(), key=lambda kv: str(kv[0])):
            names = []
            for r in runs:
                names += [p["name"] for p in r.get("phases", []) if p["name"] not in names]
            for name in names:
                t, m, _ = metric(groups, tool, dataset, task, name, b, 1)
                if t:
                    w.writerow([tool, dataset, task, b, name, f"{t['median']:.4f}", f"{t['q1']:.4f}",
                                f"{t['q3']:.4f}", f"{t['first']:.4f}",
                                f"{m['median'] / 2**30:.3f}" if m else "", t["n"]])
    with open(out / "tasks.csv", "w", newline="") as f, open(out / "outcomes.csv", "w", newline="") as g:
        w, o = csv.writer(f), csv.writer(g)
        w.writerow(["tool", "dataset", "task", "bin_size", "status", "seconds_median", "seconds_q1", "seconds_q3",
                    "seconds_first", "peak_gib_median", "peak_gib_first", "n"])
        o.writerow(["tool", "dataset", "task", "bin_size", "status", "reason"])
        for (tool, dataset, task, b), runs in sorted(groups.items(), key=lambda kv: str(kv[0])):
            t, m, status = metric(groups, tool, dataset, task, None, b, 1)
            o.writerow([tool, dataset, task, b, status, outcome(runs)[1][:300]])
            w.writerow([tool, dataset, task, b, status] + (
                [f"{t['median']:.4f}", f"{t['q1']:.4f}", f"{t['q3']:.4f}", f"{t['first']:.4f}",
                 f"{m['median'] / 2**30:.3f}" if m else "", f"{m['first'] / 2**30:.3f}" if m else "", t["n"]]
                if t else [""] * 7))


def fmt_seconds(s):
    return "—" if s is None else (f"{s * 1000:.0f} ms" if s < 1 else f"{s:.1f} s" if s < 600 else f"{s / 60:.1f} min")


def write_markdown(groups, datasets, out):
    tools = [t for t in TOOLS if any(k[0] == t for k in groups)]
    lines = ["# Benchmark summary", "",
             "Median of repeats after the first (time; peak resident memory of the whole process tree).",
             "n.a.: the tool has no equivalent; timeout/oom/fail: the run did not finish.", ""]
    for dataset in datasets:
        lines += [f"## {dataset}", "", "| Metric | " + " | ".join(LABELS[t] for t in tools) + " |",
                  "|---|" + "---|" * len(tools)]
        for title, task, phase, b, div in KEY_METRICS:
            cells = []
            for tool in tools:
                t, m, status = metric(groups, tool, dataset, task, phase, b, div)
                if t:
                    mem = f" · {m['median'] / 2**30:.1f} GiB" if m else ""
                    cells.append(f"{fmt_seconds(t['median'])}{mem}")
                else:
                    cells.append(status if status != "missing" else "")
            if any(cells):
                lines.append(f"| {title} | " + " | ".join(cells) + " |")
        lines.append("")
    (out / "summary.md").write_text("\n".join(lines))


# -- figures ---------------------------------------------------------------------------

def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
        "ytick.color": INK2, "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlecolor": INK,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False})
    return plt


def _save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", dpi=200, bbox_inches="tight")


def fig_metrics(groups, dataset, out, kind):
    """Small multiples: one panel per key metric, one bar per tool (log scale)."""
    plt = _style()
    tools = [t for t in TOOLS if any(k[0] == t for k in groups)]
    metrics = [m for m in KEY_METRICS if any((t, dataset, m[1], m[3]) in groups for t in tools)]
    if not metrics:
        return
    cols = 4
    rows = int(np.ceil(len(metrics) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(7.2, 1.55 * rows + 0.5), squeeze=False)
    for ax, (title, task, phase, b, div) in zip(axes.flat, metrics):
        ys = np.arange(len(tools))
        values = []
        for y, tool in zip(ys, tools):
            t, m, status = metric(groups, tool, dataset, task, phase, b, div)
            s = (t if kind == "time" else m)
            if s:
                v = s["median"] if kind == "time" else s["median"] / 2**30
                lo = (s["q1"] if kind == "time" else s["q1"] / 2**30)
                hi = (s["q3"] if kind == "time" else s["q3"] / 2**30)
                ax.barh(y, v, height=0.62, color=COLOURS[tool], edgecolor=SURFACE, linewidth=1)
                ax.errorbar(v, y, xerr=[[max(v - lo, 0)], [max(hi - v, 0)]], color=INK2, lw=0.8, capsize=0)
                values.append(v)
            else:
                ax.text(0.02, y, {"missing": "—"}.get(status, status), transform=ax.get_yaxis_transform(),
                        va="center", ha="left", color=INK2, fontsize=7, style="italic")
        ax.set_xscale("log")
        ax.set_yticks(ys, [LABELS[t].replace(" (ours)", "") for t in tools], fontsize=7)
        ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
        ax.set_title(title, loc="left")
        if values:
            lo, hi = min(values), max(values)
            ax.set_xlim(lo / 3, hi * 3)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    fig.supxlabel("seconds (median, IQR; log scale)" if kind == "time" else
                  "peak memory of the process tree, GiB (median, IQR; log scale)", color=INK2, fontsize=8)
    fig.tight_layout()
    _save(fig, out, f"fig_{dataset}_{kind}")
    plt.close(fig)


def fig_scaling(groups, out, kind):
    """Synthetic sweep: metric vs in-tissue 2 µm bins, one line per tool, direct labels."""
    plt = _style()
    sizes = {}
    for runs in groups.values():
        for r in runs:
            if r["dataset"].startswith("syn_") and r.get("n_bins_2um"):
                sizes[r["dataset"]] = r["n_bins_2um"]
    if not sizes:
        return
    datasets = sorted(sizes, key=sizes.get)
    tools = [t for t in TOOLS if any(k[0] == t for k in groups)]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), squeeze=False)
    for ax, (title, task, phase, b, div) in zip(axes.flat, SCALING):
        for tool in tools:
            xs, ys = [], []
            for d in datasets:
                t, m, _ = metric(groups, tool, d, task, phase, b, div)
                s = t if kind == "time" else m
                if s:
                    xs.append(sizes[d])
                    ys.append(s["median"] if kind == "time" else s["median"] / 2**30)
            if xs:
                ax.plot(xs, ys, color=COLOURS[tool], lw=2, marker="o", ms=4.5, mec=SURFACE, mew=1.5)
                ax.annotate(LABELS[tool].replace(" (ours)", ""), (xs[-1], ys[-1]), xytext=(4, 0),
                            textcoords="offset points", va="center", fontsize=6.5, color=INK2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, loc="left")
        ax.set_xlabel("in-tissue 2 µm bins")
    fig.supylabel("seconds" if kind == "time" else "peak memory (GiB)", color=INK2, fontsize=8)
    handles = [plt.Line2D([], [], color=COLOURS[t], lw=2, marker="o", ms=4.5, label=LABELS[t]) for t in tools]
    fig.legend(handles=handles, loc="upper center", ncol=len(tools), bbox_to_anchor=(0.5, 1.02), fontsize=7)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, out, f"fig_scaling_{kind}")
    plt.close(fig)


def fig_outcomes(groups, dataset, out):
    """Which tasks each tool performs on `dataset`, and how the runs ended (text in each cell)."""
    plt = _style()
    tools = [t for t in TOOLS if any(k[0] == t for k in groups)]
    columns = sorted({(k[2], k[3]) for k in groups if k[1] == dataset},
                     key=lambda c: (["image", "bins", "cells", "e2e"].index(c[0].split(".")[0]), c[0], c[1] or 0))
    if not columns:
        return
    shade = {"ok": "#dcebfa", "partial": "#f3e3c3", "n.a.": "#f1f0ed", "timeout": "#f6d4d3", "oom": "#f6d4d3",
             "fail": "#f6d4d3", "missing": SURFACE}
    fig, ax = plt.subplots(figsize=(max(7.2, 0.42 * len(columns) + 1.5), 0.42 * len(tools) + 1.6))
    for i, tool in enumerate(tools):
        for j, (task, b) in enumerate(columns):
            runs = groups.get((tool, dataset, task, b))
            status = outcome(runs)[0] if runs else "missing"
            ax.add_patch(plt.Rectangle((j, i), 0.96, 0.92, color=shade.get(status, SURFACE)))
            ax.text(j + 0.48, i + 0.46, {"ok": "✓", "missing": ""}.get(status, status), ha="center", va="center",
                    fontsize=6.5, color=INK)
    ax.set_xlim(0, len(columns))
    ax.set_ylim(len(tools), 0)
    ax.set_yticks(np.arange(len(tools)) + 0.46, [LABELS[t] for t in tools])
    ax.set_xticks(np.arange(len(columns)) + 0.48,
                  [f"{t}" + (f" {b}µm" if b else "") for t, b in columns], rotation=60, ha="right", fontsize=6.5)
    ax.grid(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(f"Task coverage and outcome · {dataset}", loc="left")
    fig.tight_layout()
    _save(fig, out, f"fig_outcomes_{dataset}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="*", help="JSONL files (default: results/*.jsonl)")
    parser.add_argument("--loupe", help="manual Loupe Browser CSV (adapters/loupe/results_template.csv format)")
    parser.add_argument("--out", default=str(BENCH / "results" / "report"))
    args = parser.parse_args()
    paths = args.results or sorted((BENCH / "results").glob("*.jsonl"))
    records = load(paths, args.loupe)
    if not records:
        raise SystemExit("No results found.")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    groups = group(records)
    datasets = sorted({r["dataset"] for r in records}, key=lambda d: (d.startswith("syn_"), d))
    write_tables(groups, out)
    write_markdown(groups, datasets, out)
    for d in datasets:
        if not d.startswith("syn_"):
            fig_metrics(groups, d, out, "time")
            fig_metrics(groups, d, out, "memory")
        fig_outcomes(groups, d, out)
    fig_scaling(groups, out, "time")
    fig_scaling(groups, out, "memory")
    print(f"{len(records)} records -> {out}")


if __name__ == "__main__":
    main()
