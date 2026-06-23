"""
Benchmark VisiumVisualizer vs Loupe Browser — puntos clave:
  1. Carga de datos Visium HD (bins 2 / 8 / 16 µm)
  2. Generación de spots V2 desde bins (2 µm y 8 µm)
  3. Segmentación celular StarDist
  4. Preprocesamiento (normalización, PCA, vecinos) — 16 µm
  5. Clustering Leiden a 3 resoluciones distintas — 16 µm

Salida: benchmarks/results/benchmark_<timestamp>.csv
"""

import sys
import os
import csv
import socket
import traceback
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

import time
import psutil

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

DATA_ROOT  = Path("/scratch/tcondegerma/data/Visium_HD_Mouse_Small_Intestine")
IMAGE_PATH = DATA_ROOT / "Visium_HD_Mouse_Small_Intestine_tissue_image.btf"
TILES_PATH = DATA_ROOT / "tiles"
SEGM_PATH  = DATA_ROOT / "StarDist"

CLUSTERING_BIN_SIZE = 16          # µm — resolución usada para clustering
LEIDEN_RESOLUTIONS  = [0.5, 1.0, 2.0]

RESULTS_DIR = REPO_ROOT / "benchmarks" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------
CSV_FIELDS = [
    "timestamp", "step", "bins_size_um", "resolution",
    "duration_s",
    "n_obs", "n_vars", "n_clusters", "n_hvgs",
    "ram_before_mb", "ram_after_mb", "peak_ram_mb",
    "hostname", "notes",
]

def _ram_mb() -> float:
    return psutil.Process().memory_info().rss / 1024 ** 2


class BenchmarkLogger:
    def __init__(self, csv_path: Path):
        self.csv_path  = csv_path
        self.hostname  = socket.gethostname()
        self._fh       = csv_path.open("w", newline="", buffering=1)
        self._writer   = csv.DictWriter(self._fh, fieldnames=CSV_FIELDS)
        self._writer.writeheader()
        self._fh.flush()
        print(f"[benchmark] Logging to: {csv_path}")

    def record(
        self, step: str, duration_s: float,
        bins_size_um=None, resolution=None,
        n_obs=None, n_vars=None, n_clusters=None, n_hvgs=None,
        ram_before_mb=None, ram_after_mb=None, peak_ram_mb=None,
        notes="",
    ):
        row = {
            "timestamp":    datetime.now().isoformat(timespec="seconds"),
            "step":         step,
            "bins_size_um": bins_size_um,
            "resolution":   resolution,
            "duration_s":   round(duration_s, 3),
            "n_obs":        n_obs,
            "n_vars":       n_vars,
            "n_clusters":   n_clusters,
            "n_hvgs":       n_hvgs,
            "ram_before_mb": round(ram_before_mb, 1)  if ram_before_mb is not None else None,
            "ram_after_mb":  round(ram_after_mb,  1)  if ram_after_mb  is not None else None,
            "peak_ram_mb":   round(peak_ram_mb,   1)  if peak_ram_mb   is not None else None,
            "hostname":     self.hostname,
            "notes":        notes,
        }
        self._writer.writerow(row)
        self._fh.flush()
        tag = f"{step}" + (f" (bins={bins_size_um}µm)" if bins_size_um else "") \
                        + (f" (r={resolution})"        if resolution    else "")
        print(f"  ✓ {tag}: {duration_s:.2f} s")

    def close(self):
        self._fh.close()


@contextmanager
def timed_step():
    """Yields (result_holder, stop_fn). Call stop_fn() to capture the time."""
    t0       = time.perf_counter()
    ram_pre  = _ram_mb()
    peak     = [ram_pre]

    class _Holder:
        elapsed = 0.0
        ram_before = ram_pre
        ram_after  = 0.0
        peak_ram   = 0.0

    h = _Holder()
    try:
        yield h
    finally:
        h.elapsed    = time.perf_counter() - t0
        h.ram_after  = _ram_mb()
        h.peak_ram   = max(h.ram_before, h.ram_after)


# ---------------------------------------------------------------------------
# Benchmark steps
# ---------------------------------------------------------------------------

def bench_load(log: BenchmarkLogger, bins_size: int):
    """Carga AnnData .h5 + posiciones .parquet para un tamaño de bin."""
    from load_data import DataSet

    step = f"load_{bins_size}um"
    print(f"\n[{step}] Cargando datos Visium HD (bins={bins_size}µm)...")

    with timed_step() as h:
        ds = DataSet(
            dir_base  = str(DATA_ROOT),
            img_file  = IMAGE_PATH.name,
            bins_size = bins_size,
        )
        _, adata, _ = ds.load()

    log.record(
        step          = step,
        duration_s    = h.elapsed,
        bins_size_um  = bins_size,
        n_obs         = adata.n_obs,
        n_vars        = adata.n_vars,
        ram_before_mb = h.ram_before,
        ram_after_mb  = h.ram_after,
        peak_ram_mb   = h.peak_ram,
        notes         = f"h5+parquet read",
    )
    return ds, adata


def bench_spots(log: BenchmarkLogger, adata, tissue_positions, bins_size: int):
    """Genera spots Visium V2 a partir de bins (solo tiene sentido para bins pequeños)."""
    from utils.bins2spots import Bins2Spots

    step = f"bins2spots_{bins_size}um"
    print(f"\n[{step}] Generando spots V2 desde bins ({bins_size}µm)...")

    with timed_step() as h:
        b2s        = Bins2Spots(adata, tissue_positions, bins_size=float(bins_size))
        spots_adata = b2s.convert(normalize=False)

    n_spots = spots_adata.n_obs if spots_adata is not None else None
    log.record(
        step          = step,
        duration_s    = h.elapsed,
        bins_size_um  = bins_size,
        n_obs         = n_spots,
        n_vars        = spots_adata.n_vars if spots_adata is not None else None,
        ram_before_mb = h.ram_before,
        ram_after_mb  = h.ram_after,
        peak_ram_mb   = h.peak_ram,
        notes         = "create_bin_squares + sjoin + aggregate",
    )
    return spots_adata


def bench_segmentation(log: BenchmarkLogger, img):
    """Segmentación celular con StarDist (carga pickle si ya existe)."""
    from image import ImageProcessor

    step = f"stardist_segmentation"
    from_cache = (SEGM_PATH / "labels_norm_v1.pickle").exists()
    print(f"\n[{step}] Segmentación StarDist "
          f"({'desde caché' if from_cache else 'inferencia completa'})...")

    with timed_step() as h:
        proc       = ImageProcessor(
            image      = img,
            tiles_path = str(TILES_PATH),
            segm_path  = str(SEGM_PATH),
        )
        _, polys   = proc.segment_image()
        cells_gdf  = proc.get_cells_gdf(polys)

    n_cells = len(cells_gdf) if cells_gdf is not None else None
    log.record(
        step          = step,
        duration_s    = h.elapsed,
        n_obs         = n_cells,
        ram_before_mb = h.ram_before,
        ram_after_mb  = h.ram_after,
        peak_ram_mb   = h.peak_ram,
        notes         = "from_cache=True" if from_cache else "full_stardist_inference",
    )
    return cells_gdf


def bench_preprocessing(log: BenchmarkLogger, sample, bins_size: int):
    """Normalización, HVGs, escala, PCA, vecinos (prerequisito de clustering)."""
    import scanpy as sc

    step = f"preprocessing_{bins_size}um"
    print(f"\n[{step}] Preprocesando (normalize→PCA→neighbors, bins={bins_size}µm)...")

    # Forzar recálculo aunque haya caché
    sample.adata_processed = None

    with timed_step() as h:
        adata_proc = sample.prepare_adata_processed(n_pcs=30, n_neighbors=15)

    n_hvgs = int(adata_proc.n_vars)
    log.record(
        step          = step,
        duration_s    = h.elapsed,
        bins_size_um  = bins_size,
        n_obs         = adata_proc.n_obs,
        n_vars        = sample.adata.n_vars,
        n_hvgs        = n_hvgs,
        ram_before_mb = h.ram_before,
        ram_after_mb  = h.ram_after,
        peak_ram_mb   = h.peak_ram,
        notes         = "filter_genes+normalize+log1p+hvg+scale+pca+neighbors",
    )


def bench_clustering(log: BenchmarkLogger, sample, resolution: float, bins_size: int):
    """Leiden clustering a una resolución. Reutiliza adata_processed ya calculado."""
    import scanpy as sc

    step = f"leiden_r{resolution}"
    print(f"\n[{step}] Clustering Leiden (r={resolution}, bins={bins_size}µm)...")

    # Aseguramos que el preprocesamiento está listo sin volver a cronometrarlo
    if sample.adata_processed is None:
        sample.prepare_adata_processed()

    with timed_step() as h:
        sc.tl.leiden(sample.adata_processed, resolution=resolution)
        # Propagar clusters a adata principal (igual que en run_clustering)
        common_idx = sample.adata.obs_names.intersection(
            sample.adata_processed.obs_names
        )
        sample.adata.obs["leiden"] = None
        sample.adata.obs.loc[common_idx, "leiden"] = (
            sample.adata_processed.obs.loc[common_idx, "leiden"]
        )
        # Genes marcadores (igual que en production)
        sc.tl.rank_genes_groups(
            sample.adata_processed, groupby="leiden", method="wilcoxon"
        )

    n_clusters = sample.adata_processed.obs["leiden"].nunique()
    log.record(
        step          = step,
        duration_s    = h.elapsed,
        bins_size_um  = bins_size,
        resolution    = resolution,
        n_obs         = sample.adata.n_obs,
        n_vars        = sample.adata.n_vars,
        n_clusters    = n_clusters,
        ram_before_mb = h.ram_before,
        ram_after_mb  = h.ram_after,
        peak_ram_mb   = h.peak_ram,
        notes         = "leiden+rank_genes_groups(wilcoxon)",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    run_id   = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = RESULTS_DIR / f"benchmark_{run_id}.csv"
    log      = BenchmarkLogger(csv_path)

    print("=" * 60)
    print(f"  VisiumHD Benchmark — {run_id}")
    print(f"  Host: {socket.gethostname()}")
    print(f"  Data: {DATA_ROOT}")
    print("=" * 60)

    # ------------------------------------------------------------------
    # 1. Carga de datos — 3 tamaños de bin
    # ------------------------------------------------------------------
    loaded = {}
    for bs in [2, 8, 16]:
        try:
            ds, adata = bench_load(log, bins_size=bs)
            loaded[bs] = (ds, adata)
        except Exception:
            print(f"  ✗ Fallo cargando {bs}µm:\n{traceback.format_exc()}")
            log.record(
                step=f"load_{bs}um", duration_s=0,
                bins_size_um=bs, notes="ERROR: " + traceback.format_exc(limit=1).strip(),
            )

    # ------------------------------------------------------------------
    # 2. Generación de spots V2 — 2 µm y 8 µm
    # ------------------------------------------------------------------
    for bs in [2, 8]:
        if bs not in loaded:
            continue
        ds, adata = loaded[bs]
        try:
            bench_spots(log, adata, ds.df_tissue_positions, bins_size=bs)
        except Exception:
            print(f"  ✗ Fallo bins2spots {bs}µm:\n{traceback.format_exc()}")
            log.record(
                step=f"bins2spots_{bs}um", duration_s=0,
                bins_size_um=bs, notes="ERROR: " + traceback.format_exc(limit=1).strip(),
            )

    # ------------------------------------------------------------------
    # 3. Segmentación celular StarDist
    # ------------------------------------------------------------------
    try:
        # Cargamos la imagen solo si no tenemos ya la segmentación en caché
        if CLUSTERING_BIN_SIZE in loaded:
            ds_ref = loaded[CLUSTERING_BIN_SIZE][0]
            img    = ds_ref.img
        else:
            # Fallback: cargar imagen directamente
            import tifffile
            img = tifffile.imread(str(IMAGE_PATH))
        bench_segmentation(log, img)
    except Exception:
        print(f"  ✗ Fallo segmentación:\n{traceback.format_exc()}")
        log.record(
            step="stardist_segmentation", duration_s=0,
            notes="ERROR: " + traceback.format_exc(limit=1).strip(),
        )

    # ------------------------------------------------------------------
    # 4. Preprocesamiento + Clustering a 3 resoluciones — bins 16 µm
    # ------------------------------------------------------------------
    if CLUSTERING_BIN_SIZE not in loaded:
        print(f"\n  ✗ Sin datos {CLUSTERING_BIN_SIZE}µm — saltando preprocesamiento y clustering.")
    else:
        from load_data import DataSet
        import geopandas as gpd
        from sample import Sample

        ds_16, adata_16 = loaded[CLUSTERING_BIN_SIZE]

        # Construir Sample (necesita un GeoDataFrame)
        # Usamos la posición de los bins directamente como geometría punto
        try:
            gdf_16 = gpd.GeoDataFrame(
                ds_16.df_tissue_positions,
                geometry=gpd.points_from_xy(
                    ds_16.df_tissue_positions["pxl_col_in_fullres"],
                    ds_16.df_tissue_positions["pxl_row_in_fullres"],
                ),
            )
            # Asegurar que el índice de gdf y adata coincidan por barcode
            if "barcode" in gdf_16.columns:
                gdf_16 = gdf_16.set_index("barcode")
            sample = Sample(gdf=gdf_16, adata=adata_16, on="index")
        except Exception:
            print(f"  ✗ Fallo al construir Sample:\n{traceback.format_exc()}")
            sample = None

        if sample is not None:
            # Preprocesamiento (una sola vez)
            try:
                bench_preprocessing(log, sample, bins_size=CLUSTERING_BIN_SIZE)
            except Exception:
                print(f"  ✗ Fallo preprocesamiento:\n{traceback.format_exc()}")
                log.record(
                    step=f"preprocessing_{CLUSTERING_BIN_SIZE}um", duration_s=0,
                    bins_size_um=CLUSTERING_BIN_SIZE,
                    notes="ERROR: " + traceback.format_exc(limit=1).strip(),
                )

            # Clustering a 3 resoluciones (adata_processed reutilizado)
            for res in LEIDEN_RESOLUTIONS:
                try:
                    bench_clustering(log, sample, resolution=res,
                                     bins_size=CLUSTERING_BIN_SIZE)
                except Exception:
                    print(f"  ✗ Fallo clustering r={res}:\n{traceback.format_exc()}")
                    log.record(
                        step=f"leiden_r{res}", duration_s=0,
                        bins_size_um=CLUSTERING_BIN_SIZE, resolution=res,
                        notes="ERROR: " + traceback.format_exc(limit=1).strip(),
                    )

    log.close()
    print(f"\n{'=' * 60}")
    print(f"  Benchmark completado. Resultados en:\n  {csv_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
