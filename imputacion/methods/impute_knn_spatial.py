r"""
impute_knn_spatial.py
=======================
METODO: kNN espacial (baseline / SUELO del benchmark).

Para cada entrada observada como CERO (sea por mascara o biologica -- el
metodo NO distingue, igual que no podria distinguir un metodo real), imputa
con la media de ese gen en los k vecinos espaciales mas cercanos que SI lo
tienen observado (>0). Si ningun vecino lo tiene, se deja en 0 (no hay
evidencia para recuperar el valor).

Es deliberadamente "ingenuo": no tiene ningun mecanismo para decidir si un
cero es dropout o biologico -- simplemente intenta rellenar todos los ceros
por igual. Por eso sirve de SUELO de referencia: cualquier metodo mas
sofisticado (SpaMean-Impute, STAGATE, GraphST, MAGIC) deberia superarlo,
sobre todo en los datasets con sparsity biologica alta (nb_s10, nb_s50),
donde rellenar TODO cero a ciegas es precisamente lo que no hay que hacer.

Uso:
    python impute_knn_spatial.py
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
K = 8
CSV_OUT = os.path.join(INDIR, "benchmark_metodos_corr_genica.csv")


# --------------------------------------------------------------------------- #
# VECINOS ESPACIALES                                                          #
# --------------------------------------------------------------------------- #
def orden_vecinos(coords):
    """Para cada spot, el orden de todos los demas por cercania (mas cerca primero)."""
    d2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    return np.argsort(d2, axis=1)


# --------------------------------------------------------------------------- #
# EL METODO                                                                   #
# --------------------------------------------------------------------------- #
def impute_knn_spatial(obs, coords, k=K):
    """
    obs    : (n_spots, n_genes) matriz OBSERVADA (con ceros, sean mascara o
             biologicos -- el metodo no conoce esa distincion, como un metodo real)
    coords : (n_spots, 2)

    Devuelve la matriz imputada (misma forma que obs).
    """
    X = obs.copy()
    orden = orden_vecinos(coords)
    ceros_i, ceros_g = np.where(obs == 0)
    for i, g in zip(ceros_i, ceros_g):
        vecinos = orden[i]
        valores = obs[vecinos, g]
        no_cero = valores[valores > 0]
        if len(no_cero) > 0:
            # usa los k primeros vecinos QUE TENGAN valor (no los k mas
            # cercanos a secas, para no diluir la estimacion con ceros)
            X[i, g] = no_cero[:k].mean()
    return X


# --------------------------------------------------------------------------- #
# METRICAS (solo sobre las posiciones enmascaradas: es lo unico verificable)   #
# --------------------------------------------------------------------------- #
def evaluar(pred_mat, gt, recalibrar=True, seed=0):
    """
    Metricas RAW y, opcionalmente, con RECALIBRACION AFIN honesta (mismo
    patron que impute_stagate.py / impute_graphst.py / impute_magic.py,
    para comparar los 5 metodos con el mismo criterio en la tabla final).
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
        coords = np.asarray(a.obsm["spatial"], dtype=float)
        gt = dict(np.load(gtf))

        pred_mat = impute_knn_spatial(obs, coords, k=K)
        m = evaluar(pred_mat, gt)
        m2 = dict(dataset=nombre, metodo=f"kNN_espacial(k={K})", **m)
        resultados.append(m2)
        print(f"{nombre:<11}{m['rmse']:>9.3f}{m['mae']:>9.3f}{m['bias']:>9.3f}"
                f"{m['pearson']:>10.4f}{m['spearman']:>11.4f}{m['r2']:>9.4f}   |"
                f"{m['rmse_recal']:>9.3f}{m['bias_recal']:>9.3f}{m['r2_recal']:>8.4f}")

    if resultados:
        guardar_csv(resultados)
        print(f"\n[save] {CSV_OUT}  (+{len(resultados)} filas)")


if __name__ == "__main__":
    main()
