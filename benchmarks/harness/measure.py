"""Parent side of a benchmark run: one fresh process per (tool, dataset, task, repeat).

The child (an adapter, run with the tool's own interpreter) times its phases
(see worker.py). This process samples the resident memory of the child's whole
process tree (Python, a uvicorn server, Chromium renderer/GPU processes, the Qt
process of napari) and enforces the timeout and the memory cap.

Memory is the sum of RSS over the process tree. Pages shared between processes
(e.g. Chromium's) are counted once per process, so browser figures are an upper
bound; with ``pss=True`` the proportional set size is sampled instead (slower,
so the sampling rate is reduced).
"""
import json
import os
import signal
import subprocess
import threading
import time

import psutil

from .worker import RESULT_PREFIX


def _tree(proc):
    try:
        return [proc, *proc.children(recursive=True)]
    except psutil.Error:
        return []


def _tree_memory(proc, pss=False):
    total = 0
    for p in _tree(proc):
        try:
            total += p.memory_full_info().pss if pss else p.memory_info().rss
        except psutil.Error:
            pass
    return total


def _kill_tree(proc):
    procs = _tree(proc)
    for p in procs:
        try:
            p.kill()
        except psutil.Error:
            pass
    psutil.wait_procs(procs, timeout=10)


class Sampler(threading.Thread):
    def __init__(self, proc, hz=20.0, pss=False, cap_bytes=None):
        super().__init__(daemon=True)
        self.proc, self.interval, self.pss, self.cap = proc, 1.0 / hz, pss, cap_bytes
        self.samples = []            # (wall time, bytes)
        self.over_cap = False
        self._halt = threading.Event()

    def run(self):
        while not self._halt.is_set():
            value = _tree_memory(self.proc, self.pss)
            self.samples.append((time.time(), value))
            if self.cap and value > self.cap:
                self.over_cap = True
                _kill_tree(self.proc)
                return
            self._halt.wait(self.interval)

    def stop(self):
        self._halt.set()
        self.join(timeout=5)

    def peak(self, t0=None, t1=None):
        values = [v for t, v in self.samples if (t0 is None or t >= t0) and (t1 is None or t <= t1)]
        return max(values) if values else None

    def before(self, t):
        values = [v for ts, v in self.samples if ts <= t]
        return values[-1] if values else (self.samples[0][1] if self.samples else None)


def run_task(cmd, *, env=None, cwd=None, timeout=1800, cap_bytes=None, hz=20.0, pss=False, log_path=None):
    """Run one adapter invocation; returns the result dict with memory per phase."""
    full_env = dict(os.environ, **(env or {}))
    log = open(log_path, "w") if log_path else subprocess.DEVNULL
    started = time.time()
    popen = subprocess.Popen(cmd, env=full_env, cwd=cwd, stdout=subprocess.PIPE, stderr=log,
                             text=True, start_new_session=True)
    proc = psutil.Process(popen.pid)
    sampler = Sampler(proc, hz=hz if not pss else min(hz, 4.0), pss=pss, cap_bytes=cap_bytes)
    sampler.start()
    status = None
    try:
        stdout, _ = popen.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        _kill_tree(proc)
        stdout, _ = popen.communicate()
        status = "timeout"
    finally:
        sampler.stop()
        if log_path:
            log.close()
    # The child's own children (a server, a browser) must not outlive it
    try:
        os.killpg(popen.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass

    result = None
    for line in (stdout or "").splitlines():
        if line.startswith(RESULT_PREFIX):
            result = json.loads(line[len(RESULT_PREFIX):])
    if sampler.over_cap:
        status = "oom"
    if result is None:
        result = dict(status=status or ("oom" if popen.returncode in (-9, 137) else "fail"),
                      reason=f"no result (exit code {popen.returncode}); see log", phases=[], info={})
    elif status:
        result["status"] = status

    result["wall_seconds"] = time.time() - started
    result["peak_tree_bytes"] = sampler.peak()
    result["memory_metric"] = "pss" if pss else "rss"
    for phase in result.get("phases", []):
        t0, t1 = phase.get("t0"), phase.get("t1")
        if t0 is None:
            continue
        base = sampler.before(t0)
        peak = sampler.peak(t0, t1)
        # Very short phases may fall between samples: use the nearest sample
        if peak is None:
            peak = sampler.before(t1)
        phase["base_tree_bytes"] = base
        phase["peak_tree_bytes"] = peak
        phase["delta_tree_bytes"] = (peak - base) if (peak is not None and base is not None) else None
    return result
