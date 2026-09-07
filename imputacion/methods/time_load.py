"""
time_load.py
=============
Cronometra SOLO la carga de los .h5ad de disco, separada de la imputacion.
Sobre el sintetico actual (400x100) sera de milisegundos; el valor es dejar
la metodologia lista para cuando se escale a Visium HD real (175k bins),
donde la carga si sera significativa.
"""
import time, csv, os
import numpy as np
import anndata as ad

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]
CSV_OUT = os.path.join(INDIR, "timings_carga.csv")
REPETICIONES = 5

filas = []
print(f"{'dataset':<11}{'t_carga_media(s)':>18}{'t_carga_std(s)':>16}{'tamano_MB':>12}")
print("-" * 57)
for nombre in DATASETS:
    h5 = os.path.join(INDIR, f"dataset_{nombre}.h5ad")
    if not os.path.exists(h5):
        print(f"  [saltado] {h5}")
        continue
    tam_mb = os.path.getsize(h5) / (1024 ** 2)
    _ = ad.read_h5ad(h5)  # warm-up (cache de disco)
    tiempos = []
    for _ in range(REPETICIONES):
        t0 = time.perf_counter()
        _ = ad.read_h5ad(h5)
        tiempos.append(time.perf_counter() - t0)
    tm, ts = float(np.mean(tiempos)), float(np.std(tiempos))
    print(f"{nombre:<11}{tm:>18.5f}{ts:>16.5f}{tam_mb:>12.3f}")
    filas.append(dict(dataset=nombre, t_carga_media_s=round(tm, 6),
                      t_carga_std_s=round(ts, 6), tamano_MB=round(tam_mb, 4),
                      n_repeticiones=REPETICIONES))

with open(CSV_OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
    w.writeheader(); w.writerows(filas)
print(f"\n[save] {CSV_OUT}")