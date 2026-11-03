"""Repeatable source-manager benchmark (synthetic data, real Bokeh protocol).

Run from the repository root:
  .venv/bin/python benchmarks/source_updates.py --output benchmarks/before.json
Each scenario/repetition uses a fresh process. Timing runs exclude tracemalloc;
allocation runs report Python/NumPy allocations separately from process RSS.
"""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import threading
import time
import tracemalloc
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def worker(kind, action, rows, memory, baseline_ref=None):
    import numpy as np
    import pandas as pd
    from bokeh.document import Document
    from bokeh.document.events import ColumnDataChangedEvent
    from bokeh.protocol import Protocol
    from hd_viewer.board.source.manager import DynamicSourceManager
    from hd_viewer.board.source.parameter import Scaler, StaticParameter

    if baseline_ref:
        # Load the original manager without changing the working tree.
        original = subprocess.check_output(
            ['git', 'show', f'{baseline_ref}:hd_viewer/board/source/manager.py'], cwd=ROOT)
        module = types.ModuleType('hd_viewer.board.source._benchmark_baseline')
        module.__package__ = 'hd_viewer.board.source'
        exec(compile(original, '<baseline manager>', 'exec'), module.__dict__)
        DynamicSourceManager = module.DynamicSourceManager

    class Sample:
        def __init__(self):
            self.fetch_lock = threading.RLock()
            self.fetches = 0
            x = np.arange(rows, dtype=np.float64)
            names = np.array([f'obs_{i}' for i in range(rows)], dtype=object)
            gene = (x % 101).copy()
            gene[::11] = np.nan
            self.frame = pd.DataFrame(dict(x=x, y=x % 307, width=2., height=2.,
                                          obs_name=names, gene=gene, other=gene * 2), index=names)
            if kind == 'cells':
                offsets = np.linspace(0., 2 * np.pi, 13)
                self.frame['x'] = [list(v + np.cos(offsets)) for v in x]
                self.frame['y'] = [list(v + np.sin(offsets)) for v in x % 307]

        def fetch_data(self, colnames=None, bbox=None, **kwargs):
            self.fetches += 1
            start, end = (0, rows) if bbox is None else bbox[:2]
            return self.frame.iloc[int(start):int(end)][
                ['x', 'y', 'width', 'height', 'obs_name', *(colnames or [])]].copy()

        def is_continous(self, colname):
            return True

        def get_range(self, colname):
            return 0., 200.

    sample = Sample()
    manager = DynamicSourceManager(sample)
    manager.add_parameter(Scaler(0., 1., 'gene', attr='alpha', not_sequenced_val=0.))
    manager.add_parameter(StaticParameter('#123456', not_sequenced_val='#888888'))
    doc = Document()
    doc.add_root(manager.source)
    initial_rows = rows * 2 // 3
    manager.bbox = [0, initial_rows, 0, rows]
    if action != 'initial':
        manager.fetch_data()
    if action == 'stream':
        manager.bbox = [initial_rows, rows, 0, rows]
    to_remove = manager.indexes[::4]
    sample.fetches = 0
    events = []

    def capture(event):
        # The server serializes source changes as they occur. Freeze the data
        # now: otherwise an empty-source event can serialize the final source
        # again after a subsequent stream, double-counting its payload.
        if isinstance(event, ColumnDataChangedEvent) and event.data is None:
            current = getattr(event.model, event.attr)
            columns = current if event.cols is None else event.cols
            event.data = {name: current[name].copy() for name in columns if name in current}
        events.append(event)

    doc.on_change(capture)
    gc.collect()
    rss_before = int(Path('/proc/self/statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
    if memory:
        tracemalloc.start()
    started = time.perf_counter()
    if action in ('initial', 'stream'):
        manager.fetch_data()
    elif action == 'prune':
        manager.remove_indcs(to_remove)
    elif action == 'style':
        manager.add_parameter(StaticParameter('#654321', not_sequenced_val='#aaaaaa'))
    elif action == 'scale':
        manager.update_param('alpha', new_maxv=.7)
    elif action == 'gene':
        manager.add_parameter(Scaler(0., 1., 'other', attr='alpha', not_sequenced_val=0.))
    elapsed_ms = (time.perf_counter() - started) * 1e3
    allocations_peak = tracemalloc.get_traced_memory()[1] if memory else None
    if memory:
        tracemalloc.stop()
    rss_after = int(Path('/proc/self/statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    started = time.perf_counter()
    message = Protocol().create('PATCH-DOC', events)
    json_bytes = len(message.content_json.encode())
    buffer_bytes = sum(len(buffer.to_bytes()) for buffer in message.buffers)
    serialize_ms = (time.perf_counter() - started) * 1e3
    data = {key: np.asarray(value).tolist() if key not in ('x', 'y') or kind != 'cells'
            else [list(v) for v in value] for key, value in manager.source.data.items()}
    digest = hashlib.sha256(json.dumps(data, sort_keys=True, allow_nan=True).encode()).hexdigest()
    return dict(ms=elapsed_ms, serialize_ms=serialize_ms, allocations_peak_bytes=allocations_peak,
                rss_before_bytes=rss_before, rss_after_bytes=rss_after, peak_rss_bytes=peak_rss,
                patch_bytes=json_bytes + buffer_bytes, patch_json_bytes=json_bytes,
                patch_buffer_bytes=buffer_bytes, events=len(events), fetches=sample.fetches,
                loaded_rows=len(manager.indexes), digest=digest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--rows', type=int, default=60000)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--worker', nargs=2, metavar=('KIND', 'ACTION'))
    parser.add_argument('--memory', action='store_true')
    parser.add_argument('--baseline-ref', help='Read the original manager from this git revision')
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(*args.worker, args.rows, args.memory, args.baseline_ref)))
        return
    import bokeh
    import numpy
    import pandas
    results = dict(rows=args.rows, repeats=args.repeats, python=sys.version,
                   baseline_ref=args.baseline_ref,
                   versions=dict(bokeh=bokeh.__version__, numpy=numpy.__version__, pandas=pandas.__version__),
                   scenarios={})
    for kind in ('bins', 'cells'):
        for action in ('initial', 'stream', 'prune', 'style', 'scale', 'gene'):
            command = [sys.executable, str(Path(__file__).resolve()), '--worker', kind, action,
                       '--rows', str(args.rows)]
            if args.baseline_ref:
                command.extend(['--baseline-ref', args.baseline_ref])
            runs = [json.loads(subprocess.check_output(command, text=True, cwd=ROOT))
                    for _ in range(args.repeats)]
            allocation_run = json.loads(subprocess.check_output([*command, '--memory'], text=True, cwd=ROOT))
            result = {key: statistics.median(run[key] for run in runs)
                      for key in runs[0] if isinstance(runs[0][key], (float, int))}
            result.update(allocations_peak_bytes=allocation_run['allocations_peak_bytes'],
                          digest=runs[0]['digest'], runs=runs, allocation_run=allocation_run)
            results['scenarios'][f'{kind}/{action}'] = result
            print(f"{kind}/{action}: {result['ms']:.1f} ms, "
                  f"allocations {result['allocations_peak_bytes'] / 2**20:.1f} MiB, "
                  f"patch {result['patch_bytes'] / 2**20:.2f} MiB", flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
