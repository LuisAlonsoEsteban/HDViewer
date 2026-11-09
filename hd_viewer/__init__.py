import os
import tempfile

# scanpy compiles its kernels with numba and caches them. Next to a read-only
# install numba finds no cache location and the import fails ("cannot cache
# function ... no locator available"). Give it a writable folder before numba
# is imported. Do not use NUMBA_DISABLE_JIT instead: scanpy's kernels (e.g. in
# normalize_total) fail when run as plain Python.
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), f"visiumvisualizer-numba-{os.getuid()}"))
