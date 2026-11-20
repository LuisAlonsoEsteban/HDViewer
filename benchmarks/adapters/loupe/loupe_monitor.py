"""Time and memory of Loupe Browser tasks, measured by hand (Windows or macOS).

Loupe Browser is a closed desktop application with no scripting interface, so
each task is started and stopped by the operator (PROTOCOL.md). This script
samples the resident memory of every Loupe Browser process (the main process
and its renderer/GPU helpers) ten times a second and records, per task:
seconds between the two key presses and the peak memory in between.

    pip install psutil
    python loupe_monitor.py --dataset mouse_brain --out loupe_results.csv

At each prompt type the task row (e.g. "bins.render view:overview 8") and press
Enter as you start the action; press Enter again as soon as the view has
finished drawing. Repeat each task five times (the first counts as "first").
Type "q" to finish. Rows are appended to --out in the format of
results_template.csv, which `benchmarks/report.py --loupe` reads.
"""
import argparse
import csv
import os
import threading
import time

import psutil

NAMES = ("loupe",)       # process name fragments (case-insensitive): "Loupe Browser", "Loupe Browser Helper"


def loupe_processes():
    found = []
    for p in psutil.process_iter(["name"]):
        name = (p.info.get("name") or "").lower()
        if any(n in name for n in NAMES):
            found.append(p)
    return found


def tree_rss():
    total = 0
    for p in loupe_processes():
        try:
            total += p.memory_info().rss
        except psutil.Error:
            pass
    return total


class Sampler(threading.Thread):
    def __init__(self, hz=10.0):
        super().__init__(daemon=True)
        self.interval, self.samples, self.halt = 1.0 / hz, [], threading.Event()

    def run(self):
        while not self.halt.is_set():
            self.samples.append((time.perf_counter(), tree_rss()))
            self.halt.wait(self.interval)

    def peak(self, t0, t1):
        values = [v for t, v in self.samples if t0 <= t <= t1]
        return max(values) if values else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", default="loupe_results.csv")
    args = parser.parse_args()
    if not loupe_processes():
        print("Loupe Browser is not running: start it first (the idle memory is recorded as 'baseline').")
    sampler = Sampler()
    sampler.start()
    new = not os.path.exists(args.out)
    counts = {}
    with open(args.out, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["dataset", "task", "bin_size", "phase", "repeat", "seconds", "peak_rss_mb", "status", "notes"])
        print(f"Idle Loupe Browser memory: {tree_rss() / 2**20:.0f} MiB")
        while True:
            line = input("\ntask phase [bin_size]  (e.g. 'bins.render view:overview 8', q to quit): ").strip()
            if line.lower() in ("q", "quit", "exit"):
                break
            parts = line.split()
            if len(parts) < 2:
                print("Give at least a task and a phase.")
                continue
            task, phase = parts[0], parts[1]
            bin_size = parts[2] if len(parts) > 2 else ""
            input("Press Enter as you START the action...")
            t0 = time.perf_counter()
            input("Press Enter when the view has FINISHED drawing...")
            t1 = time.perf_counter()
            key = (task, phase, bin_size)
            repeat = counts.get(key, 0)
            counts[key] = repeat + 1
            notes = input("Notes (Enter for none; 'fail' if it did not finish): ").strip()
            status = ("fail" if notes.lower().startswith("fail") else
                      "n.a." if notes.lower().startswith("n.a") else "ok")
            peak = sampler.peak(t0, t1) / 2**20
            w.writerow([args.dataset, task, bin_size, phase, repeat, f"{t1 - t0:.3f}", f"{peak:.0f}", status, notes])
            f.flush()
            print(f"  {task} {phase} #{repeat}: {t1 - t0:.2f} s, peak {peak:.0f} MiB")
    sampler.halt.set()


if __name__ == "__main__":
    main()
