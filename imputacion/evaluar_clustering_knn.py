"""
evaluar_clustering_knn.py
==========================
Extiende cluster_domains.py: mide si aplicar kNN espacial (impute_knn_spatial)
sobre la matriz observada mejora la recuperacion de dominios via Leiden,
comparado con el baseline SIN imputar, Y ADEMAS cuanto acierta el metodo al
reconstruir los valores concretos que se le escondieron.

Desde que gen_datasets3.py genera con mask_frac=0.10, cada dataset lleva un
held-out real: un 10% de los genes de cada spot esta puesto a 0 en X y sus
valores verdaderos viven en ground_truth_<nombre>.npz. Eso permite cerrar el
experimento que se planteo originalmente, con las dos mitades a la vez:

  1) IMPUTACION (¿recupera los VALORES?) -- RMSE / MAE / bias / Pearson /
     Spearman / R2 sobre las posiciones enmascaradas, exactamente el mismo
     calculo de impute_knn_spatial.evaluar que produce benchmark_metodos.csv.
  2) DOMINIOS (¿recupera la ESTRUCTURA?) -- ARI/NMI/AMI/HOMO de Leiden contra
     obs["dominio"], sobre esa misma matriz imputada.

Cruzando 1) y 2) se responde la pregunta real: ¿un metodo que recupera bien los
valores (el kNN llegaba a R2~0.47 en poisson) recupera tambien mejor los
dominios que uno que no?

Tres filas de clustering por dataset, para poder leer el margen:
    sin_heldout  layers["ground_truth_full"] -- counts sin enmascarar: TECHO,
                 lo maximo que podria conseguir una imputacion perfecta
    sin_imputar  a.X                         -- dato degradado: SUELO
    knn_spatial  kNN espacial sobre a.X      -- lo que consigue el metodo

OJO con la naturaleza de los ceros: en X conviven los ceros del held-out (los
unicos verificables) y los ceros BIOLOGICOS de la binomial negativa (~10% en
nb_s10, ~0% en nb_s00). impute_knn_spatial rellena TODOS los ceros por igual,
sin distinguirlos -- es deliberado (es el suelo del benchmark), pero explica
que en nb_s10 pueda empeorar: reescribe tambien ceros que eran correctos.

Corre sobre nb_s00_dom_libre (sparsity biologica ~0%: la degradacion es casi
solo el held-out, margen limpio) y nb_s10_dom_libre (sparsity biologica ~10%
encima del held-out: escenario exigente).

Reutiliza sin modificar:
    - construir_grafo / elegir_resolucion / evaluar / cargar   de cluster_domains.py
    - impute_knn_spatial / evaluar (metricas de imputacion)    de impute_knn_spatial.py

Se recalculan TODAS las filas en esta misma ejecucion, con el MISMO seed y
MISMO barrido de resolucion -- evita depender de metricas_clustering.csv (que
quedaria obsoleto si se regeneran los datasets) y garantiza que "antes" y
"despues" se comparan de forma honesta.

Salidas:
    datasets_gen3/benchmark_clustering_dominios.csv   (ARI/NMI/AMI/HOMO)
    datasets_gen3/benchmark_imputacion_dominios.csv   (RMSE/Pearson/R2 en held-out)

Uso:
    python evaluar_clustering_knn.py
"""

import os
import warnings
import numpy as np
import pandas as pd

from cluster_domains import (OUTDIR, CFG, cargar, construir_grafo,
                            elegir_resolucion, evaluar, datasets_disponibles)

# impute_knn_spatial vive en methods/, no aqui: se anade esa carpeta al path
# para no depender de que PYTHONPATH la incluya
import sys
from rutas import METHODS, ruta_ground_truth
sys.path.insert(0, METHODS)
from impute_knn_spatial import (
    impute_knn_spatial, evaluar as evaluar_imputacion, K as K_KNN,
)

# Por defecto TODAS las variantes "libre" que haya en datasets_gen3/ (una por
# cada base de gen_datasets2: gp_crudo, gaussian, poisson, nb_s00, nb_s10,
# nb_s50), descubiertas con datasets_disponibles() -- antes eran solo las dos
# NB porque eran las unicas que gen_datasets3 generaba. Las "atada" quedan
# fuera del default: "libre" es el modo que build_domain_offsets documenta como
# recomendado (aisla la senal de dominio de la de modulo); se anaden con
# --datasets. Este script no tiene --datasets: para restringir la lista, editala aqui.
DATASETS = datasets_disponibles(patron="_dom_libre")
CSV_OUT = os.path.join(OUTDIR, "benchmark_clustering_dominios.csv")
CSV_IMP = os.path.join(OUTDIR, "benchmark_imputacion_dominios.csv")


def correr_pipeline(a, dom_real, n_target, cfg, etiqueta_metodo, nombre_dataset):
    """Corre construir_grafo + elegir_resolucion + evaluar (de cluster_domains.py,
    sin modificar) sobre el AnnData dado (a.X = la matriz a evaluar)."""
    b = construir_grafo(a, cfg)
    res, leiden, nclust, barrido = elegir_resolucion(b, n_target, cfg)
    m = evaluar(dom_real, leiden)
    m.update(
        dataset=nombre_dataset,
        metodo=etiqueta_metodo,
        n_dominios=n_target,
        sparsity_pct=round(100 * float((np.asarray(a.X) == 0).mean()), 2),
        resolution=res,
    )
    barrido_str = " ".join(f"r{r}->{n}c" for r, n in barrido)
    print(f"[{nombre_dataset:>18} | {etiqueta_metodo:>11}]  "
          f"sparsity={m['sparsity_pct']:>5.2f}%  res={res}  n_clust={m['n_clusters']:>2}  "
          f"ARI={m['ARI']:.3f}  NMI={m['NMI']:.3f}  AMI={m['AMI']:.3f}  HOMO={m['HOMO']:.3f}")
    print(f"{'':>38}barrido: {barrido_str}")
    return m


def cargar_held_out(nombre):
    """Lee ground_truth_<nombre>.npz (posiciones enmascaradas + valor real).
    Devuelve None si el dataset se genero sin held-out (mask_frac=0): en ese
    caso el .npz existe pero con arrays vacios y no hay nada que evaluar."""
    ruta = ruta_ground_truth(nombre, OUTDIR)
    if not os.path.exists(ruta):
        return None
    gt = dict(np.load(ruta))
    if len(gt["true_value"]) == 0:
        return None
    return gt


def main():
    warnings.filterwarnings("ignore")
    # --- verificacion de sanidad: held-out real presente en la generacion ---
    for nombre in DATASETS:
        a_chequeo = cargar(nombre)
        mf = a_chequeo.uns["generation_params"].get("mask_frac")
        n_mask = int(np.asarray(a_chequeo.layers["mask"]).sum())
        print(f"[chequeo] {nombre}: mask_frac={mf}  entradas enmascaradas={n_mask}")
        if n_mask == 0:
            print("          AVISO: sin held-out -> solo se podra medir ARI, no RMSE/R2. "
                  "Regenera con mask_frac>0 en gen_datasets3.CFG3.")

    filas, filas_imp = [], []

    for nombre in DATASETS:
        a = cargar(nombre)
        dom_real = a.obs["dominio"].astype(int).values
        coords = np.asarray(a.obsm["spatial"], dtype=float)
        n_target = int(len(np.unique(dom_real)))
        gt = cargar_held_out(nombre)

        # --- techo: los counts SIN enmascarar (lo maximo recuperable) ---
        a_techo = a.copy()
        a_techo.X = np.asarray(a.layers["ground_truth_full"], dtype=np.float32)
        filas.append(correr_pipeline(a_techo, dom_real, n_target, CFG,
                                     etiqueta_metodo="sin_heldout",
                                     nombre_dataset=nombre))

        # --- suelo: SIN imputar, la matriz observada tal cual ---
        filas.append(correr_pipeline(a, dom_real, n_target, CFG,
                                     etiqueta_metodo="sin_imputar",
                                     nombre_dataset=nombre))

        # --- kNN espacial: imputa, luego el MISMO pipeline de clustering ---
        obs = np.asarray(a.X, dtype=float)
        pred = impute_knn_spatial(obs, coords, k=K_KNN)

        # 1) ¿recupera los VALORES? -- solo sobre las posiciones enmascaradas
        if gt is not None:
            m_imp = evaluar_imputacion(pred, gt)
            m_imp = dict(dataset=nombre, metodo=f"kNN_espacial(k={K_KNN})",
                         n_heldout=int(len(gt["true_value"])), **m_imp)
            filas_imp.append(m_imp)
            print(f"[{nombre:>18} | {'imputacion':>11}]  "
                  f"n={m_imp['n_heldout']}  RMSE={m_imp['rmse']:.3f}  "
                  f"MAE={m_imp['mae']:.3f}  bias={m_imp['bias']:+.3f}  "
                  f"pearson={m_imp['pearson']:.4f}  R2={m_imp['r2']:.4f}  "
                  f"| R2_recal={m_imp['r2_recal']:.4f}")

        # 2) ¿recupera la ESTRUCTURA? -- clustering sobre la matriz imputada
        a_imputado = a.copy()
        a_imputado.X = pred.astype(np.float32)
        filas.append(correr_pipeline(a_imputado, dom_real, n_target, CFG,
                                     etiqueta_metodo="knn_spatial",
                                     nombre_dataset=nombre))

    # ------------------------------------------------------------------ #
    # TABLA 1: recuperacion de DOMINIOS                                   #
    # ------------------------------------------------------------------ #
    cols = ["dataset", "metodo", "sparsity_pct", "n_dominios", "n_clusters",
            "ARI", "NMI", "AMI", "HOMO", "resolution"]
    df = pd.DataFrame(filas)[cols].round(4)
    df.to_csv(CSV_OUT, index=False)
    print(f"\nTabla clustering -> {CSV_OUT}")
    print(df.to_string(index=False))

    piv = df.pivot(index="dataset", columns="metodo", values="ARI")
    if {"sin_imputar", "knn_spatial"}.issubset(piv.columns):
        piv["delta_ARI"] = (piv["knn_spatial"] - piv["sin_imputar"]).round(4)
    if {"sin_imputar", "sin_heldout"}.issubset(piv.columns):
        # margen disponible: cuanto ARI perdio el held-out y podria devolverse
        piv["margen"] = (piv["sin_heldout"] - piv["sin_imputar"]).round(4)
        if "delta_ARI" in piv.columns:
            # fraccion del margen que el metodo recupera de verdad
            piv["frac_recuperada"] = (piv["delta_ARI"] / piv["margen"].replace(0, np.nan)).round(4)
    print("\nARI por metodo (techo sin_heldout, suelo sin_imputar, metodo knn_spatial):")
    print(piv.to_string())

    # ------------------------------------------------------------------ #
    # TABLA 2: recuperacion de VALORES en el held-out                     #
    # ------------------------------------------------------------------ #
    if filas_imp:
        df_imp = pd.DataFrame(filas_imp).round(4)
        df_imp.to_csv(CSV_IMP, index=False)
        print(f"\nTabla imputacion (held-out) -> {CSV_IMP}")
        cols_imp = ["dataset", "metodo", "n_heldout", "rmse", "mae", "bias",
                    "pearson", "spearman", "r2", "rmse_recal", "r2_recal"]
        print(df_imp[[c for c in cols_imp if c in df_imp.columns]].to_string(index=False))

        # --- el cruce: ¿recuperar valores implica recuperar dominios? ---
        print("\nCruce valores vs dominios (por dataset):")
        for nombre in DATASETS:
            fila_imp = df_imp[df_imp["dataset"] == nombre]
            if fila_imp.empty or nombre not in piv.index:
                continue
            r2 = float(fila_imp["r2"].iloc[0])
            print(f"  {nombre:>18}:  R2_heldout={r2:>7.4f}   "
                  f"ARI suelo={piv.loc[nombre, 'sin_imputar']:.4f} -> "
                  f"knn={piv.loc[nombre, 'knn_spatial']:.4f} "
                  f"(techo {piv.loc[nombre, 'sin_heldout']:.4f})")
    else:
        print("\n[!] Ningun dataset tenia held-out: no hay metricas de imputacion. "
              "Regenera con mask_frac=0.10 en gen_datasets3.CFG3.")


if __name__ == "__main__":
    main()
