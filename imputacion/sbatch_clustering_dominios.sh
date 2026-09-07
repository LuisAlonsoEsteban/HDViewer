"""
evaluar_clustering_metodos.py  (version con --metodo)
========================================================
Uso:
    python evaluar_clustering_metodos.py --metodo knn
    python evaluar_clustering_metodos.py --metodo spamean
    python evaluar_clustering_metodos.py --metodo stagate
    python evaluar_clustering_metodos.py --metodo magic
    python evaluar_clustering_metodos.py               (sin --metodo: los 4, como antes)
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats

from cluster_domains import OUTDIR, CFG, cargar, construir_grafo, elegir_resolucion, evaluar
from rutas import METHODS
sys.path.insert(0, METHODS)

from time_methods import _run_knn, _run_spamean, _run_stagate, _run_magic

DATASETS = ["nb_s00_dom_libre", "nb_s10_dom_libre"]

# clave corta (para --metodo y para el nombre de fichero) -> (nombre bonito, funcion)
METODOS_DISPONIBLES = {
    "knn":     ("kNN_espacial(k=8)", _run_knn),
    "spamean": ("SpaMean",           _run_spamean),
    "stagate": ("STAGATE",           _run_stagate),
    "magic":   ("MAGIC",             _run_magic),
}


def calc_metricas_imputacion(pred_mat, gt):
    pred = pred_mat[gt["spot_idx"], gt["gene_idx"]]
    true = gt["true_value"]
    err = pred - true
    rmse = float(np.sqrt(np.mean(err ** 2)))
    mae = float(np.mean(np.abs(err)))
    bias = float(np.mean(err))
    pearson = float(stats.pearsonr(pred, true)[0]) if np.std(pred) > 0 else float("nan")
    spearman = float(stats.spearmanr(pred, true)[0]) if np.std(pred) > 0 else float("nan")
    ss_res = float(np.sum(err ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    rng = np.random.default_rng(0)
    idx = rng.permutation(len(true))
    mitad = len(idx) // 2
    cal_idx, test_idx = idx[:mitad], idx[mitad:]
    a, b = np.polyfit(pred[cal_idx], true[cal_idx], deg=1)
    pred_cal = a * pred[test_idx] + b
    err_cal = pred_cal - true[test_idx]
    r2_recal = 1.0 - np.sum(err_cal**2) / np.sum((true[test_idx] - true[test_idx].mean())**2)

    return dict(n_heldout=len(true), rmse=rmse, mae=mae, bias=bias,
               pearson=pearson, spearman=spearman, r2=r2, r2_recal=float(r2_recal))


def correr_clustering(a, dom_real, n_target, cfg, etiqueta, nombre_dataset):
    b = construir_grafo(a, cfg)
    res, leiden, nclust, barrido = elegir_resolucion(b, n_target, cfg)
    m = evaluar(dom_real, leiden)
    m.update(dataset=nombre_dataset, metodo=etiqueta, n_dominios=n_target,
             sparsity_pct=round(100 * float((np.asarray(a.X) == 0).mean()), 2),
             resolution=res)
    print(f"[{nombre_dataset:>18} | {etiqueta:>18}]  sparsity={m['sparsity_pct']:>5.2f}%  "
          f"res={res}  n_clust={m['n_clusters']:>2}  ARI={m['ARI']:.3f}  HOMO={m['HOMO']:.3f}")
    return m


def main(metodo_clave):
    warnings.filterwarnings("ignore")

    # --- seleccion de metodo(s): uno solo, o los 4 si no se especifica ---
    if metodo_clave is None:
        metodos_activos = METODOS_DISPONIBLES
        sufijo = ""  # ficheros como antes, sin sufijo (comportamiento por defecto)
    else:
        if metodo_clave not in METODOS_DISPONIBLES:
            validos = ", ".join(METODOS_DISPONIBLES.keys())
            raise ValueError(f"--metodo '{metodo_clave}' no reconocido. Validos: {validos}")
        metodos_activos = {metodo_clave: METODOS_DISPONIBLES[metodo_clave]}
        sufijo = f"_{metodo_clave}"  # ej. "_knn", "_stagate"

    csv_cluster = os.path.join(OUTDIR, f"benchmark_clustering_dominios{sufijo}.csv")
    csv_imput = os.path.join(OUTDIR, f"benchmark_imputacion_dominios{sufijo}.csv")

    filas_cluster, filas_imput = [], []

    for nombre in DATASETS:
        a = cargar(nombre)
        dom_real = a.obs["dominio"].astype(int).values
        n_target = int(len(np.unique(dom_real)))
        gt = np.load(os.path.join(OUTDIR, f"ground_truth_{nombre}.npz"))

        # techo y suelo se recalculan SIEMPRE, para que cada CSV por metodo
        # sea autocontenido (no dependa de otro fichero para interpretar frac_recuperada)
        a_full = a.copy()
        a_full.X = np.asarray(a.layers["ground_truth_full"], dtype=np.float32)
        m_techo = correr_clustering(a_full, dom_real, n_target, CFG, "sin_heldout", nombre)
        filas_cluster.append(m_techo)

        m_suelo = correr_clustering(a, dom_real, n_target, CFG, "sin_imputar", nombre)
        filas_cluster.append(m_suelo)

        for clave, (nombre_bonito, run_fn) in metodos_activos.items():
            pred = run_fn(a)

            m_imput = calc_metricas_imputacion(pred, gt)
            m_imput.update(dataset=nombre, metodo=nombre_bonito)
            filas_imput.append(m_imput)

            a_imputado = a.copy()
            a_imputado.X = pred.astype(np.float32)
            m_clust = correr_clustering(a_imputado, dom_real, n_target, CFG,
                                        nombre_bonito, nombre)
            m_clust["frac_recuperada"] = (m_clust["ARI"] / m_techo["ARI"]
                                          if m_techo["ARI"] > 0 else float("nan"))
            filas_cluster.append(m_clust)

    cols_c = ["dataset", "metodo", "sparsity_pct", "n_dominios", "n_clusters", "ARI", "NMI", "AMI", "HOMO", "resolution", "frac_recuperada"]
    df_c = pd.DataFrame(filas_cluster)
    for c in cols_c:
        if c not in df_c.columns:
            df_c[c] = np.nan
    df_c = df_c[cols_c].round(4)
    df_c.to_csv(csv_cluster, index=False)

    cols_i = ["dataset", "metodo", "n_heldout", "rmse", "mae", "bias", "pearson", "spearman", "r2", "r2_recal"]
    df_i = pd.DataFrame(filas_imput)[cols_i].round(4)
    df_i.to_csv(csv_imput, index=False)

    print(f"\nTabla clustering -> {csv_cluster}")
    print(df_c.to_string(index=False))
    print(f"\nTabla imputacion (held-out) -> {csv_imput}")
    print(df_i.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metodo", type=str, default=None,
                        choices=list(METODOS_DISPONIBLES.keys()),
                        help="Ejecutar un solo metodo (knn/spamean/stagate/magic). "
                             "Si se omite, se ejecutan los 4.")
    args = parser.parse_args()
    main(args.metodo)