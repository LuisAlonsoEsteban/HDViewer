"""
time_methods.py
=================
Cronometra la EJECUCION de cada metodo de imputacion de forma justa y
reproducible, SEPARADA de la evaluacion de metricas (que se hace en los
impute_*.py). Mide por separado:
  - t_carga    : leer el .h5ad de disco
  - t_imputacion : SOLO la llamada al metodo (impute_X), sin metricas ni I/O

Hace 1 warm-up descartado (absorbe imports en frio, cache, JIT) + N
repeticiones cronometradas; reporta media y desviacion estandar.

Guarda en un CSV APARTE (timings_metodos.csv), nunca mezclado con
benchmark_metodos.csv (metricas).

Uso:
    python time_methods.py --metodo knn      [--repeticiones 5]
    python time_methods.py --metodo spamean
    python time_methods.py --metodo stagate
    python time_methods.py --metodo graphst
    python time_methods.py --metodo magic
"""

import argparse
import time
import csv
import os
import numpy as np
import anndata as ad

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]
CSV_TIMINGS = os.path.join(INDIR, "timings_metodos.csv")


# --------------------------------------------------------------------------- #
# ADAPTADORES: cada metodo con su firma real, devolviendo SOLO la imputacion  #
# --------------------------------------------------------------------------- #
def _run_knn(a):
    from impute_knn_spatial import impute_knn_spatial, K
    obs = np.asarray(a.X, dtype=float)
    coords = np.asarray(a.obsm["spatial"], dtype=float)
    return impute_knn_spatial(obs, coords, k=K)


def _run_spamean(a):
    from impute_spamean import impute_spamean, K, DELTA
    obs = np.asarray(a.X, dtype=float)
    coords = np.asarray(a.obsm["spatial"], dtype=float)
    pred, _ = impute_spamean(obs, coords, k=K, delta=DELTA)
    return pred


def _run_stagate(a):
    from impute_stagate import impute_stagate
    pred, _ = impute_stagate(a)
    return pred

def _run_graphst(a):
    from impute_graphst import impute_graphst
    pred, _, _ = impute_graphst(a, verbose=False)
    return pred


def _run_magic(a):
    from impute_magic import impute_magic
    obs = np.asarray(a.X, dtype=float)
    return impute_magic(obs)


ADAPTADORES = {
    "knn": ("kNN_espacial", _run_knn),
    "spamean": ("SpaMean-Impute", _run_spamean),
    "stagate": ("STAGATE_pyG", _run_stagate),
    "graphst": ("GraphST", _run_graphst),
    "magic": ("MAGIC", _run_magic),
}


# --------------------------------------------------------------------------- #
# CRONOMETRAJE                                                                 #
# --------------------------------------------------------------------------- #
def cronometrar(metodo_key, repeticiones):
    metodo_nombre, run_fn = ADAPTADORES[metodo_key]
    filas = []

    print(f"{'dataset':<11}{'t_carga(s)':>14}{'t_imput_media(s)':>18}"
          f"{'t_imput_std(s)':>16}{'n_rep':>7}")
    print("-" * 66)

    for nombre in DATASETS:
        h5 = os.path.join(INDIR, f"dataset_{nombre}.h5ad")
        if not os.path.exists(h5):
            print(f"  [saltado] no existe {h5}")
            continue

        # --- t_carga: leer el .h5ad (medido una vez, es determinista) ---
        t0 = time.perf_counter()
        a = ad.read_h5ad(h5)
        t_carga = time.perf_counter() - t0

        # --- warm-up descartado (absorbe imports en frio / cache / JIT) ---
        _ = run_fn(a)

        # --- N repeticiones cronometradas de SOLO la imputacion ---
        tiempos = []
        for _ in range(repeticiones):
            t0 = time.perf_counter()
            _ = run_fn(a)
            tiempos.append(time.perf_counter() - t0)

        t_media = float(np.mean(tiempos))
        t_std = float(np.std(tiempos))

        print(f"{nombre:<11}{t_carga:>14.4f}{t_media:>18.4f}{t_std:>16.4f}"
              f"{repeticiones:>7}")

        filas.append(dict(
            metodo=metodo_nombre, dataset=nombre,
            t_carga_s=round(t_carga, 6),
            t_imputacion_media_s=round(t_media, 6),
            t_imputacion_std_s=round(t_std, 6),
            n_repeticiones=repeticiones,
            tiempos_crudos_s=";".join(f"{t:.6f}" for t in tiempos),
        ))

    # --- guardar en CSV APARTE, con append para acumular todos los metodos ---
    existe = os.path.exists(CSV_TIMINGS)
    with open(CSV_TIMINGS, "a" if existe else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        if not existe:
            w.writeheader()
        w.writerows(filas)
    print(f"\n[save] {CSV_TIMINGS}  (+{len(filas)} filas para {metodo_nombre})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--metodo", required=True, choices=list(ADAPTADORES.keys()))
    ap.add_argument("--repeticiones", type=int, default=5)
    args = ap.parse_args()
    cronometrar(args.metodo, args.repeticiones)