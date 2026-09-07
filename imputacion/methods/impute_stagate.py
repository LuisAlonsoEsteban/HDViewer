r"""
impute_stagate.py
====================
METODO: STAGATE_pyG (Dong & Zhang, 2022, Nat Commun) -- autoencoder de
atencion sobre grafo espacial. A diferencia de MAGIC (grafo de expresion) o
kNN/SpaMean (heuristicas fijas), aqui la reconstruccion en espacio de genes
sale de un modelo ENTRENADO (grafo espacial + autoencoder).

DECISION METODOLOGICA (importante, dejar constancia en la memoria):
STAGATE normalmente se alimenta con datos normalizados y log-transformados
(normalize_total + log1p), como en su tutorial estandar. Aqui, en cambio,
se entrena con LOS COUNTS CRUDOS, sin normalizar, para que la reconstruccion
quede en la MISMA ESCALA que el ground truth y sea directamente comparable
con kNN/SpaMean/MAGIC (que tampoco normalizan). Es una desviacion deliberada
del tutorial estandar -- valida para este benchmark sintetico porque el
proceso generativo no tiene efectos de batch ni de profundidad de secuenciacion
que "normalizar" tuviera que corregir. Queda como ablacion pendiente probar
la version log-normalizada si se quiere comparar.

rad_cutoff: se usa 3.0 porque coincide con el length_scale=3.0 del generador
(gen_datasets.py) -- da un vecindario espacial acorde a la escala de
correlacion real de los datos.

Requiere: STAGATE_pyG instalado (ya clonado/instalado en methods/STAGATE_pyG)

Uso:
    python impute_stagate.py
"""

import numpy as np
import anndata as ad
import scanpy as sc
from scipy import stats
import os
import csv

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]
RAD_CUTOFF = 3.0    # NO usar el 150 del tutorial: ese valor es para
                    # coordenadas Visium reales en micrometros (~100um entre
                    # spots). Nuestro grid sintetico tiene espaciado=1, asi
                    # que 3.0 (coherente con length_scale=3 del generador)
                    # es el equivalente correcto.
N_EPOCHS = 3000
PREPROCESAMIENTO = "crudo"   # "crudo" | "normalizado" (estilo tutorial)
TARGET_SUM = 1e4            # solo si PREPROCESAMIENTO="normalizado"
CSV_OUT = os.path.join(INDIR, "benchmark_metodos.csv")


# --------------------------------------------------------------------------- #
# EL METODO                                                                    #
# --------------------------------------------------------------------------- #
def impute_stagate(adata, rad_cutoff=RAD_CUTOFF, seed=0, n_epochs=N_EPOCHS,
                   preprocesamiento=PREPROCESAMIENTO):
    """
    adata : AnnData con .X = counts CRUDOS y .obsm['spatial'].

    preprocesamiento="crudo": entrena directo sobre counts (decision original,
        evita tener que invertir ninguna transformacion, pero puede converger
        peor numericamente -- ver nota de bias sistematico en el chat).

    preprocesamiento="normalizado": replica el tutorial oficial de STAGATE
        (normalize_total + log1p ANTES de construir el grafo y entrenar).
        Tras reconstruir, se DESHACE la transformacion para volver a la
        escala de counts crudos y poder comparar contra el ground truth:
            1. X_norm_recon = expm1(ReX)                    -- deshace log1p
            2. X_raw_recon  = X_norm_recon * (lib_size_i / target_sum)
                                                              -- deshace normalize_total
        usando el tamano de libreria (suma de counts por spot) ORIGINAL,
        capturado antes de normalizar.

    Devuelve la matriz reconstruida en ESCALA DE COUNTS CRUDOS en ambos
    casos, para que sea comparable de forma identica al resto de metodos.
    """
    import STAGATE_pyG

    adata = adata.copy()
    adata.X = adata.X.astype(np.float32)

    lib_size = None
    if preprocesamiento == "normalizado":
        lib_size = adata.X.sum(axis=1).copy()          # tamano de libreria ORIGINAL
        lib_size = np.asarray(lib_size).ravel()
        lib_size = np.clip(lib_size, 1.0, None)          # evita division por 0
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        sc.pp.log1p(adata)
    elif preprocesamiento != "crudo":
        raise ValueError(f"preprocesamiento desconocido: {preprocesamiento!r}")

    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff)
    STAGATE_pyG.Stats_Spatial_Net(adata)

    adata = STAGATE_pyG.train_STAGATE(adata, random_seed=seed, n_epochs=n_epochs, save_reconstrction=True)   # sic: typo real

    candidatos_layers = ["STAGATE_ReX", "ReX"]
    candidatos_obsm = ["STAGATE_ReX", "ReX"]

    rex = None
    for k in candidatos_layers:
        if k in adata.layers:
            rex = np.asarray(adata.layers[k]); break
    if rex is None:
        for k in candidatos_obsm:
            if k in adata.obsm:
                rex = np.asarray(adata.obsm[k]); break

    if rex is None:
        print("\n[AVISO] No se encontro la reconstruccion en las claves esperadas.")
        print(f"  adata.layers disponibles: {list(adata.layers.keys())}")
        print(f"  adata.obsm disponibles:   {list(adata.obsm.keys())}")
        print("  Ajusta 'candidatos_layers'/'candidatos_obsm' con el nombre correcto.\n")
        raise KeyError("Reconstruccion de STAGATE no encontrada; ver aviso arriba.")

    if preprocesamiento == "normalizado":
        # deshacer log1p y luego deshacer normalize_total, por spot
        rex = np.expm1(rex)
        rex = rex * (lib_size[:, None] / TARGET_SUM)

    return rex, adata


# --------------------------------------------------------------------------- #
# METRICAS (identicas al resto de scripts de metodos)                         #
# --------------------------------------------------------------------------- #
def evaluar(pred_mat, gt, recalibrar=True, seed=0):
    """
    Metricas RAW y, opcionalmente, con RECALIBRACION AFIN honesta:
    y = a*pred + b, ajustada con la MITAD de las entradas enmascaradas y
    evaluada en la OTRA MITAD (holdout) -- para no hacer trampa ajustando
    y evaluando sobre los mismos datos.

    Util cuando el bias es sistematico (ver nota en impute_stagate.py):
    si la correlacion (Pearson/Spearman) ya es razonable pero el bias es
    grande, es un problema de ESCALA, no de que el modelo no capture la
    señal -- la recalibracion lo demuestra directamente.
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

        # ajusta y = a*pred + b SOLO con la mitad de calibracion
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
        gt = dict(np.load(gtf))

        pred_mat, _ = impute_stagate(a)
        m = evaluar(pred_mat, gt)
        m2 = dict(dataset=nombre, metodo=f"STAGATE_pyG_{PREPROCESAMIENTO}", **m)
        resultados.append(m2)
        print(f"{nombre:<11}{m['rmse']:>9.3f}{m['mae']:>9.3f}{m['bias']:>9.3f}"
              f"{m['pearson']:>10.4f}{m['spearman']:>11.4f}{m['r2']:>9.4f}   |"
              f"{m['rmse_recal']:>9.3f}{m['bias_recal']:>9.3f}{m['r2_recal']:>8.4f}")

    if resultados:
        guardar_csv(resultados)
        print(f"\n[save] {CSV_OUT}  (+{len(resultados)} filas)")
        print("\nColumnas '_r' = tras RECALIBRACION AFIN honesta (y=a*pred+b, ajustada")
        print("con la mitad de las entradas de mascara, evaluada en la otra mitad).")
        print("Si bias_r ~ 0 y r2_r mejora mucho respecto a r2, el modelo SI capta la")
        print("señal espacial pero devuelve la escala equivocada (problema de")
        print("convergencia numerica con counts crudos, no de que 'no aprenda nada').")


if __name__ == "__main__":
    main()
