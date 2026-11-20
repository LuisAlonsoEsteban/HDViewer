"""In-process side of a benchmark run (imported by every adapter, in every tool env).

Standard library only, so it works in any interpreter. An adapter receives the
task, dataset and options as JSON on the command line, records named *phases*,
and prints one ``@@RESULT {json}`` line on stdout. The parent process
(``harness/measure.py``) samples the RSS of the whole process tree and assigns
peak memory to each phase using the wall-clock bounds recorded here.
"""
import argparse
import contextlib
import gc
import json
import os
import platform
import sys
import time
import traceback

RESULT_PREFIX = "@@RESULT "


class NotApplicable(Exception):
    """The tool has no equivalent of this task (recorded as status 'n.a.')."""


class Recorder:
    def __init__(self):
        self.phases = []
        self.info = {}
        self._tracemalloc = os.environ.get("BENCH_TRACEMALLOC") == "1"

    @contextlib.contextmanager
    def phase(self, name, **extra):
        """Time a block. `extra` and anything added to the yielded dict are stored."""
        gc.collect()
        record = dict(name=name, **extra)
        if self._tracemalloc:
            import tracemalloc
            tracemalloc.start()
        record["t0"] = time.time()
        started = time.perf_counter()
        try:
            yield record
        finally:
            record["seconds"] = time.perf_counter() - started
            record["t1"] = time.time()
            if self._tracemalloc:
                import tracemalloc
                record["py_alloc_peak_bytes"] = tracemalloc.get_traced_memory()[1]
                tracemalloc.stop()
            self.phases.append(record)

    def add(self, name, seconds, **extra):
        """Record a phase timed elsewhere (e.g. inside a browser)."""
        now = time.time()
        self.phases.append(dict(name=name, seconds=seconds, t0=now - seconds, t1=now, **extra))

    def note(self, **values):
        self.info.update(values)


def _versions(modules):
    out = {}
    for name in modules:
        try:
            from importlib.metadata import version
            out[name] = version(name)
        except Exception:
            pass
    return out


def main(tasks, versions=()):
    """Entry point of an adapter: ``python adapters/x.py --task ID --spec JSON``.

    `tasks` maps task id -> function(spec, rec). A function raises NotApplicable
    when the tool cannot do the task.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", required=True)
    parser.add_argument("--spec", required=True, help="JSON: dataset + task options")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    spec = json.loads(args.spec)
    rec = Recorder()
    out = dict(task=args.task, status="ok", python=platform.python_version(), pid=os.getpid(),
               versions=_versions(versions))
    try:
        func = tasks.get(args.task)
        if func is None:
            raise NotApplicable("task not implemented by this tool")
        func(spec, rec)
    except NotApplicable as exc:
        out.update(status="n.a.", reason=str(exc))
    except MemoryError:
        out.update(status="oom", reason="MemoryError")
    except Exception as exc:  # recorded, never dropped
        out.update(status="fail", reason=f"{type(exc).__name__}: {exc}",
                   traceback=traceback.format_exc()[-4000:])
    out.update(phases=rec.phases, info=rec.info)
    sys.stdout.flush()
    print(RESULT_PREFIX + json.dumps(out, default=str), flush=True)


def wait_until(predicate, timeout=600.0, interval=0.01):
    """Poll `predicate` until true; raises TimeoutError."""
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise TimeoutError("condition not reached")
        time.sleep(interval)
