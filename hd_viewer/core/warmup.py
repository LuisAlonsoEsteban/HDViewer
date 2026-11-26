"""Analysis runtime set-up at server start.

scanpy runs its neighbour search (PyNNDescent) and t-SNE single-threaded unless
told otherwise, and numba compiles PyNNDescent/UMAP kernels on first use (~10 s).
`warm_up` sets the thread count and runs a tiny clustering pass in the background,
so the first visitor's analysis does not pay the compilation.
"""
import os
import threading
import time

from .. import settings


def analysis_threads():
    if settings.ANALYSIS_THREADS > 0:
        return settings.ANALYSIS_THREADS
    try:
        return len(os.sched_getaffinity(0))          # respects cgroup/SLURM CPU limits
    except AttributeError:
        return os.cpu_count() or 1


def configure_scanpy():
    import scanpy as sc
    sc.settings.n_jobs = analysis_threads()
    return sc


def _warm_up():
    started = time.monotonic()
    try:
        import anndata as ad
        import numpy as np

        sc = configure_scanpy()
        rng = np.random.default_rng(0)
        adata = ad.AnnData(rng.poisson(1.0, size=(300, 40)).astype(np.float32))
        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
        print(f"[startup] Analysis kernels ready in {time.monotonic() - started:.1f} s")
    except Exception as exc:                          # never blocks the server: only slower first runs
        print(f"[startup] Analysis warm-up skipped: {exc}")


_started = None
_lock = threading.Lock()


def warm_up():
    """Start the warm-up in a daemon thread (once per process); returns it."""
    global _started
    with _lock:
        if _started is None:
            _started = threading.Thread(target=_warm_up, daemon=True, name="analysis-warmup")
            _started.start()
        return _started
