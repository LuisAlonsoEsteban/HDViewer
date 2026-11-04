"""Benchmark viewport lookup and complete fetches using the real sample backend.

Run with --baseline-ref to load the original backend read-only from Git.
Each scenario/repetition runs in a fresh process. No browser or network is used.
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
import time
import tracemalloc
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SCENARIOS = {
    'bins_native_256': ('bins', 256, 1, 0., 80.),
    'bins_native_768': ('bins', 768, 1, 0., 80.),
    'cells_native_768': ('cells', 768, 1, 0., 80.),
    'bins_aggregate4_768': ('bins', 768, 4, 0., 80.),
    'bins_rotated_aggregate4_768': ('bins', 768, 4, 27., 80.),
    'bins_aggregate4_wide_768': ('bins', 768, 4, 0., 320.),
}


def rss():
    return int(Path('/proc/self/statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')


def measure(operation, allocations=False):
    if allocations:
        tracemalloc.start()
    started = time.perf_counter()
    result = operation()
    elapsed = (time.perf_counter() - started) * 1000
    peak = tracemalloc.get_traced_memory()[1] if allocations else None
    if allocations:
        tracemalloc.stop()
    return result, elapsed, peak


def worker(name, baseline_ref, allocations):
    # Match the app's initialization order; startup loads UI/backend packages.
    import hd_viewer.startup
    import anndata as ad
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import shapely
    from hd_viewer.startup.data.image.transform import SpatialTransform
    from hd_viewer.startup.sample import Bins, Cells

    if baseline_ref:
        for filename, alias in (('_sample', '_viewport_baseline_sample'),
                                ('bins', '_viewport_baseline_bins'),
                                ('cells', '_viewport_baseline_cells')):
            source = subprocess.check_output([
                'git', 'show', f'{baseline_ref}:hd_viewer/startup/sample/{filename}.py'], cwd=ROOT).decode()
            source = source.replace('from ._sample import Sample',
                                    'from ._viewport_baseline_sample import Sample')
            module = types.ModuleType(f'hd_viewer.startup.sample.{alias}')
            module.__package__ = 'hd_viewer.startup.sample'
            sys.modules[module.__name__] = module
            exec(compile(source, f'<baseline {filename}>', 'exec'), module.__dict__)
        Bins = sys.modules['hd_viewer.startup.sample._viewport_baseline_bins'].Bins
        Cells = sys.modules['hd_viewer.startup.sample._viewport_baseline_cells'].Cells

    kind, side, factor, angle_degrees, window = SCENARIOS[name]
    row, col = np.indices((side, side))
    row, col = row.ravel(), col.ravel()
    theta = np.deg2rad(angle_degrees)
    matrix = np.array([[2 * np.cos(theta), -2 * np.sin(theta), side],
                       [2 * np.sin(theta), 2 * np.cos(theta), side], [0., 0., 1.]])
    identity = tuple(np.eye(3).ravel())
    transform = SpatialTransform(identity, identity, tuple(matrix.ravel()), side * 5,
                                 1., 2., f'benchmark-{angle_degrees}')
    xy = transform.array_points(np.column_stack((col, row)))
    if kind == 'cells':
        geometries = shapely.box(xy[:, 0] - .7, xy[:, 1] - .6, xy[:, 0] + .7, xy[:, 1] + .6)
    else:
        geometries = transform.rectangles(xy)
    names = pd.Index([f'b{i}' for i in range(len(row))])
    gdf = gpd.GeoDataFrame(dict(index=names, array_col=col, array_row=row),
                           geometry=geometries, index=names)
    counts = (np.arange(len(row), dtype=np.float32) % 101).reshape(-1, 1)
    adata = ad.AnnData(counts, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=['gene']))
    sample = (Bins(gdf, adata, bin_size=2, pixels_per_bin=2, spatial_transform=transform)
              if kind == 'bins' else Cells(gdf, adata))
    level = f'agg:{factor}' if factor > 1 else None
    # Prepare common caches equally; measure the added query-index cost below.
    sample.get_coordinates(gdf)
    sample.get_obs_names()
    sample._gene_vector('gene')
    if factor > 1:
        sample.aggregation_grid(factor)
    centre = transform.array_points([[side * .5, side * .5]])[0]
    boxes = [[centre[0] + i * 1.37 - window/2, centre[0] + i * 1.37 + window/2,
              centre[1] + i * .91 - window/2, centre[1] + i * .91 + window/2] for i in range(9)]

    def lookup(bbox):
        if factor == 1:
            return sample.crop_geom(*bbox)
        if hasattr(sample, '_crop_aggregation_grid'):
            cropped, _ = sample._crop_aggregation_grid(factor, bbox)
            return cropped
        # Original lookup kernel from Bins._fetch_aggregated_data, excluding
        # the subsequent groupby/value preparation measured by fetch below.
        grid = sample.aggregation_grid(factor)
        left, right, bottom, top = bbox
        angle = sample.spatial_transform.angle
        half_x = (abs(np.cos(angle))*sample.native_width*factor + abs(np.sin(angle))*sample.native_height*factor)/2
        half_y = (abs(np.sin(angle))*sample.native_width*factor + abs(np.cos(angle))*sample.native_height*factor)/2
        mask = ((grid.x + half_x >= left) & (grid.x - half_x <= right)
                & (grid.y + half_y >= bottom) & (grid.y - half_y <= top))
        cropped = sample.gdf.loc[mask]
        grid.loc[mask].copy()
        return cropped

    def fetch(bbox):
        if factor > 1:
            sample._multiscale_cache.clear()  # measure uncached viewports
        return sample.fetch_data(bbox=bbox, colnames=['gene'], resolution_level=level)

    gc.collect()
    rss_before = rss()
    _, cold_ms, cold_peak = measure(lambda: lookup(boxes[0]), allocations)
    gc.collect()
    rss_indexed = rss()
    lookup_times, fetch_times, lookup_peaks, fetch_peaks = [], [], [], []
    digest = hashlib.sha256()
    selected_rows = []
    for bbox in boxes[1:]:
        selected, elapsed, peak = measure(lambda: lookup(bbox), allocations)
        lookup_times.append(elapsed)
        lookup_peaks.append(peak or 0)
        selected_rows.append(len(selected))
        data, elapsed, peak = measure(lambda: fetch(bbox), allocations)
        fetch_times.append(elapsed)
        fetch_peaks.append(peak or 0)
        digest.update(data.to_json(orient='split', double_precision=15).encode())
        del selected, data
    return dict(rows=len(gdf), selected_rows=selected_rows,
                cold_lookup_ms=cold_ms, lookup_ms=statistics.median(lookup_times),
                fetch_ms=statistics.median(fetch_times),
                cold_peak_bytes=cold_peak, lookup_peak_bytes=max(lookup_peaks),
                fetch_peak_bytes=max(fetch_peaks), rss_before_bytes=rss_before,
                rss_indexed_bytes=rss_indexed, rss_after_bytes=rss(),
                process_peak_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                digest=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-ref')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--worker', choices=SCENARIOS)
    parser.add_argument('--allocations', action='store_true')
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.baseline_ref, args.allocations)))
        return
    results = dict(baseline_ref=args.baseline_ref, repeats=args.repeats, scenarios={})
    for name in SCENARIOS:
        command = [sys.executable, str(Path(__file__).resolve()), '--worker', name]
        if args.baseline_ref:
            command.extend(['--baseline-ref', args.baseline_ref])
        runs = [json.loads(subprocess.check_output(command, cwd=ROOT, text=True)) for _ in range(args.repeats)]
        allocation_run = json.loads(subprocess.check_output([*command, '--allocations'], cwd=ROOT, text=True))
        result = {key: statistics.median(run[key] for run in runs)
                  for key in runs[0] if isinstance(runs[0][key], (int, float))}
        for key in ('cold_peak_bytes', 'lookup_peak_bytes', 'fetch_peak_bytes'):
            result[key] = allocation_run[key]
        result.update(digest=runs[0]['digest'], selected_rows=runs[0]['selected_rows'],
                      runs=runs, allocation_run=allocation_run)
        results['scenarios'][name] = result
        print(f"{name}: cold {result['cold_lookup_ms']:.1f} ms, lookup {result['lookup_ms']:.2f} ms, "
              f"fetch {result['fetch_ms']:.2f} ms, lookup allocation {result['lookup_peak_bytes']/2**20:.2f} MiB", flush=True)
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
