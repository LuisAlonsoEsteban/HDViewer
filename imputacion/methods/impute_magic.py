r"""
impute_magic.py
=================
METODO: MAGIC (van Dijk et al., 2018, Cell) -- via el paquete magic-impute.

A diferencia de kNN espacial y SpaMean-Impute, MAGIC NO usa coordenadas
espaciales en absoluto: construye un grafo de vecinos en el espacio de
EXPRESION (que spots tienen perfiles de genes parecidos, no que spots estan
cerca fisicamente) y difunde valores sobre ese grafo con un operador de
Markov. Sirve de referencia "sin informacion espacial" -- si STAGATE/GraphST
no le sacan ventaja clara, la senal espacial no estaria aportando gran cosa.

MAGIC tampoco tiene ningun mecanismo para distinguir dropout de cero
biologico: como SoftImpute/SimpleImpute/scVI/gimVI/tangram en el paper de
Hafiz et al. (2026), rellena TODOS los ceros por igual.

Requiere: pip install magic-impute  (ya instalado en tfm-impute)

Uso (desde imputacion/ o imputacion/methods/, la ruta es absoluta):
    python impute_magic.py
"""

import numpy as np
import anndata as ad
from scipy import stats
import os
import csv

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]
CSV_OUT = os.path.join(INDIR, "benchmark_metodos.csv")


# --------------------------------------------------------------------------- #
# EL METODO                                                                    #
# --------------------------------------------------------------------------- #
def impute_magic(obs, seed=0):
    """
    obs : (n_spots, n_genes) matriz OBSERVADA (con ceros, sin distinguir origen)
    MAGIC no usa coordenadas -- solo la matriz de expresion.
    Devuelve la matriz imputada (misma forma).
    """
    import magic
    magic_op = magic.MAGIC(random_state=seed, verbose=False)
    X_imputado = magic_op.fit_transform(obs)
    return np.asarray(X_imputado)


# --------------------------------------------------------------------------- #
# METRICAS (identicas a impute_knn_spatial.py / impute_spamean.py)             #
# --------------------------------------------------------------------------- #
def evaluar(pred_mat, gt, recalibrar=True, seed=0):
    """
    Metricas RAW y, opcionalmente, con RECALIBRACION AFIN honesta:
    y = a*pred + b, ajustada con la MITAD de las entradas enmascaradas y
    evaluada en la OTRA MITAD (holdout) -- mismo patron que impute_stagate.py
    e impute_graphst.py, para comparar los 4 metodos con el mismo criterio.
    """
    pred_full = pred_mat[gt["spot_idx"], gt["gene_idx"]]
    true_full = gt["true_value"]

    def calc(pred, true):
        err = pred - true
        rmse = float(np.sqrt(np.mean(err ** 2)))
        mae = float(np.mean(np.abs(err)))
        bias = float(np.mean(err))
        pearson = float(stats.pearsonr(pred, true)[0]) if np.std(pred) > 0 else float("nan")
        spearman = float(stats.spearmanr(pred, true)[0]) if np.std(pred) > 0 else float("nan")
        ss_res = float(np.sum(err ** 2))
        ss_tot = float(np.sum((true - true.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return dict(rmse=rmse, mae=mae, bias=bias, pearson=pearson, spearman=spearman, r2=r2)

    resultado = calc(pred_full, true_full)

    if recalibrar:
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(true_full))
        mitad = len(idx) // 2
        cal_idx, test_idx = idx[:mitad], idx[mitad:]

        a, b = np.polyfit(pred_full[cal_idx], true_full[cal_idx], deg=1)
        pred_cal = a * pred_full[test_idx] + b

        resultado_cal = calc(pred_cal, true_full[test_idx])
        resultado_cal = {f"{k}_recal": v for k, v in resultado_cal.items()}
        resultado.update(resultado_cal)
        resultado["recal_a"] = float(a)
        resultado["recal_b"] = float(b)

    return resultado

def guardar_csv(filas, path=CSV_OUT):
    existe = os.path.exists(path)
    campos = list(filas[0].keys())
    with open(path, "a" if existe else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        if not existe:
            w.writeheader()
        w.writerows(filas)


# --------------------------------------------------------------------------- #
# MAIN                                                                         #
# --------------------------------------------------------------------------- #
def main():
    resultados = []
    print(f"{'dataset':<11}{'RMSE':>9}{'MAE':>9}{'bias':>9}{'pearson':>10}"
      f"{'spearman':>11}{'r2':>9}   |{'RMSE_r':>9}{'bias_r':>9}{'r2_r':>8}")
    print("-" * 100)
    for nombre in DATASETS:
        h5 = os.path.join(INDIR, f"dataset_{nombre}.h5ad")
        gtf = os.path.join(INDIR, f"ground_truth_{nombre}.npz")
        if not os.path.exists(h5):
            print(f"  [saltado] no existe {h5}")
            continue
        a = ad.read_h5ad(h5)
        obs = np.asarray(a.X, dtype=float)
        gt = dict(np.load(gtf))

        pred_mat = impute_magic(obs)
        m = evaluar(pred_mat, gt)
        m2 = dict(dataset=nombre, metodo="MAGIC", **m)
        resultados.append(m2)
        print(f"{nombre:<11}{m['rmse']:>9.3f}{m['mae']:>9.3f}{m['bias']:>9.3f}"
                f"{m['pearson']:>10.4f}{m['spearman']:>11.4f}{m['r2']:>9.4f}   |"
                f"{m['rmse_recal']:>9.3f}{m['bias_recal']:>9.3f}{m['r2_recal']:>8.4f}")

    if resultados:
        guardar_csv(resultados)
        print(f"\n[save] {CSV_OUT}  (+{len(resultados)} filas)")


if __name__ == "__main__":
    main()