r"""
impute_spamean.py
===================
METODO: SpaMean-Impute (replica de Hafiz, Azim & Shatabda, 2026,
"Spatial information matters: are traditional imputation methods effective
for spatial transcriptomics data?", Briefings in Bioinformatics).

A diferencia del kNN espacial (impute_knn_spatial.py), que rellena TODO cero
por igual, SpaMean-Impute tiene un mecanismo EXPLICITO para decidir si un
cero es dropout (recuperable) o biologico (dejar como esta):

    1. Para cada spot i, se buscan sus k vecinos espaciales mas cercanos
       N_i^(k) (excluyendo el propio spot).
    2. Para cada gen j con X_ij = 0, se calcula el "zero_ratio":
           zero_ratio_ij = (# vecinos de i con el gen j tambien en 0) / k
       Es decir: si CASI TODOS los vecinos tambien tienen ese gen en 0,
       el ratio es alto -> parece zona donde el gen genuinamente no se
       expresa (biologico). Si POCOS vecinos lo tienen en 0, el ratio es
       bajo -> el gen SI se expresa alrededor, esto es sospechoso de dropout.
    3. Umbral de decision, el "dropout threshold" delta:
           zero_ratio_ij <= delta  ->  DROPOUT VALIDO -> se imputa
           zero_ratio_ij >  delta  ->  BIOLOGICO       -> se deja en 0
    4. La imputacion es la media de los valores NO CERO de ese gen entre
       los k vecinos (si ninguno es no-cero, no se puede imputar: se deja).

Parametros (ver "Average zero sparsity..." y sensibilidad del paper):
    k     : num. vecinos espaciales. El paper reporta buen comportamiento
            entre 3 y 11. Aqui se usa 8 (igual que el resto del proyecto).
    delta : umbral de dropout. Rango util en el paper: 4%-40%. Con delta bajo
            el metodo es CONSERVADOR (imputa pocos ceros, mas falsos negativos
            de dropout); con delta alto imputa mas (mas riesgo de tocar
            ceros biologicos). Aqui se usa 0.4 (el ejemplo del paper).

NOTA (para la mejora futura que se comento): este metodo usa una regla FIJA
(umbral simple sobre zero_ratio) y una imputacion muy simple (media). Es el
punto de partida natural para mejorarlo -- p.ej. aprender el umbral, ponderar
por distancia en vez de voto uniforme, o sustituir la media por un modelo.

Uso:
    python impute_spamean.py
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
DELTA = 0.4
CSV_OUT = os.path.join(INDIR, "benchmark_metodos.csv")


# --------------------------------------------------------------------------- #
# VECINOS ESPACIALES (k-mas-cercanos fijos, no todo el orden)                  #
# --------------------------------------------------------------------------- #
def knn_indices(coords, k):
    """Para cada spot, los indices de sus k vecinos mas cercanos (sin el propio)."""
    d2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    return np.argsort(d2, axis=1)[:, :k]


# --------------------------------------------------------------------------- #
# EL METODO (Fig. 2 del paper)                                                #
# --------------------------------------------------------------------------- #
def impute_spamean(obs, coords, k=K, delta=DELTA):
    """
    obs    : (n_spots, n_genes) matriz OBSERVADA
    coords : (n_spots, 2)

    Devuelve (imputado, zero_ratio_registrado) -- zero_ratio_registrado es un
    dict {(i,g): ratio} solo para las entradas que eran cero, util para
    inspeccionar despues cuantas se consideraron dropout vs biologicas.
    """
    X = obs.copy()
    vecinos = knn_indices(coords, k)
    ceros_i, ceros_g = np.where(obs == 0)

    n_dropout_detectado = 0
    n_biologico_detectado = 0

    for i, g in zip(ceros_i, ceros_g):
        vec = vecinos[i]
        valores_vecinos = obs[vec, g]

        zero_ratio = float(np.mean(valores_vecinos == 0))   # ecuacion del paper

        if zero_ratio <= delta:
            # dropout valido -> imputar con la media de los vecinos NO cero
            no_cero = valores_vecinos[valores_vecinos > 0]
            if len(no_cero) > 0:
                X[i, g] = no_cero.mean()
                n_dropout_detectado += 1
            # si todos los vecinos son 0 pese a pasar el umbral (raro,
            # posible con k pequeno), no hay con que imputar: se deja
        else:
            # considerado biologicamente no expresado -> NO TOCAR
            n_biologico_detectado += 1

    return X, dict(n_dropout_detectado=n_dropout_detectado,
                   n_biologico_detectado=n_biologico_detectado)


# --------------------------------------------------------------------------- #
# METRICAS (igual que en impute_knn_spatial.py, para comparar en la misma tabla)
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
      f"{'spearman':>11}{'r2':>9}   |{'RMSE_r':>9}{'bias_r':>9}{'r2_r':>8}"
      f"   {'n_drop':>8}{'n_bio':>8}")
    print("-" * 110)
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

        pred_mat, info = impute_spamean(obs, coords, k=K, delta=DELTA)
        m = evaluar(pred_mat, gt)
        m2 = dict(dataset=nombre, metodo=f"SpaMean-Impute(k={K},delta={DELTA})", **m)
        resultados.append(m2)
        print(f"{nombre:<11}{m['rmse']:>9.3f}{m['mae']:>9.3f}{m['bias']:>9.3f}"
                f"{m['pearson']:>10.4f}{m['spearman']:>11.4f}{m['r2']:>9.4f}   |"
                f"{m['rmse_recal']:>9.3f}{m['bias_recal']:>9.3f}{m['r2_recal']:>8.4f}"
                f"   {info['n_dropout_detectado']:>8}{info['n_biologico_detectado']:>8}")

    if resultados:
        guardar_csv(resultados)
        print(f"\n[save] {CSV_OUT}  (+{len(resultados)} filas)")


if __name__ == "__main__":
    main()
