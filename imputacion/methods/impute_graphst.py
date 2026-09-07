r"""
impute_graphst.py
====================
METODO: GraphST (Long et al., 2023, Nat Commun) -- aprendizaje contrastivo
auto-supervisado sobre grafo espacial.

AVISO METODOLOGICO IMPORTANTE:
El output principal de GraphST es un EMBEDDING de baja dimension
(adata.obsm['emb']), pensado para clustering/dominios espaciales, NO
necesariamente una matriz reconstruida en el espacio original de genes
comparable directamente contra el ground truth de counts.

Este script localiza una reconstruccion en espacio de genes tras entrenar
e IMPRIME EXPLICITAMENTE que clave uso -- imprescindible para confiar en el
resultado: un R2~0 puede significar "la reconstruccion es real pero sin
senal" (resultado legitimo) o "se engancho una clave que no es una
reconstruccion real" (bug). Sin saber la clave, no se puede distinguir.

Uso:
    python impute_graphst.py             -> los 4 datasets de gen_datasets2
                                            (synthetic_out/datasets/)
    python impute_graphst.py --gen3      -> los 4 datasets con dominios de
                                            gen_datasets3 (synthetic_out/datasets_gen3/)
    python impute_graphst.py --barrido   -> barrido n_neighbors x dim_output x alpha

Las rutas SIEMPRE salen de rutas.py, nunca relativas al CWD: da igual desde
donde se lance el script.
"""

import numpy as np
import anndata as ad
from scipy import stats
import os
import csv
import torch
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "GraphST"))

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR, GEN3, ruta_dataset, ruta_ground_truth
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]
CSV_OUT = os.path.join(INDIR, "benchmark_metodos.csv")

# --gen3: los datasets con DOMINIOS espaciales de gen_datasets3.py. Desde que
# se generan con mask_frac=0.10 llevan held-out real (ground_truth_*.npz ya no
# esta vacio), asi que se pueden evaluar con EXACTAMENTE el mismo criterio que
# los de arriba -- y luego cruzar ese R2 con el ARI de dominios.
DATASETS_GEN3 = ["nb_s00_dom_libre", "nb_s00_dom_atada",
                 "nb_s10_dom_libre", "nb_s10_dom_atada"]
CSV_OUT_GEN3 = os.path.join(GEN3, "benchmark_metodos_gen3.csv")


# --------------------------------------------------------------------------- #
# EL METODO                                                                    #
# --------------------------------------------------------------------------- #
def impute_graphst(adata, seed=42, n_neighbors=3, dim_output=64, alpha=10, beta=1,
                   verbose=True):
    """
    alpha: peso de la perdida de RECONSTRUCCION (reconstruir la expresion
           original a partir del embedding) en la funcion de perdida total.
    beta : peso de la perdida CONTRASTIVA (acercar/alejar pares positivos y
           negativos en el embedding) en la funcion de perdida total.
    Ojo: no se ha podido verificar contra el codigo fuente instalado si
    "mas alpha" implica mas peso relativo real o si hay normalizacion interna
    -- tratar como hiperparametros a barrer empiricamente, no a razonar
    a priori sobre su efecto exacto.

    Devuelve (matriz_reconstruida_o_None, adata_entrenado, clave_usada_o_None).
    """
    from GraphST import GraphST as GraphST_mod
    from GraphST.preprocess import construct_interaction

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    adata = adata.copy()
    adata.X = adata.X.astype(np.float32)

    construct_interaction(adata, n_neighbors=n_neighbors)

    model = GraphST_mod.GraphST(
        adata, device=device, random_seed=seed,
        dim_output=dim_output, alpha=alpha, beta=beta,
    )
    adata = model.train()

    n_genes = adata.shape[1]
    # OJO: 'emb' se ha quitado de los candidatos. Coincidia en FORMA con
    # n_genes solo por casualidad (tu dataset tiene 100 genes), pero es un
    # espacio latente aprendido -- NO una reconstruccion de expresion.
    # Prueba empirica: con 'emb' el bias era ~ -mean(counts) en los 4
    # datasets (embedding centrado en 0, sin relacion de escala con los
    # counts reales) y la recalibracion afin no arreglaba nada (r2_recal~0)
    # -- confirma que no hay correspondencia recuperable, no solo una mala
    # calibracion de escala como si pasaba con STAGATE.
    candidatos = ["recon", "X_rec", "reconstruction", "denoised", "emb"]

    if verbose:
        print(f"  [diagnostico] adata.obsm disponibles: "
              f"{[(k, v.shape) for k, v in adata.obsm.items()]}")
        print(f"  [diagnostico] adata.layers disponibles: {list(adata.layers.keys())}")
        print(f"  [diagnostico] n_genes esperado (para que la clave sirva "
              f"como reconstruccion): {n_genes}")
        if "emb" in adata.obsm and "emb" not in candidatos:
            print(f"  [diagnostico] 'emb' existe (shape={adata.obsm['emb'].shape}) ")

    for k in candidatos:
        if k in adata.obsm and adata.obsm[k].shape[1] == n_genes:
            if verbose:
                print(f"  [diagnostico] CLAVE USADA: adata.obsm['{k}']  "
                      f"shape={adata.obsm[k].shape}")
            return np.asarray(adata.obsm[k]), adata, k
    for k in candidatos:
        if k in adata.layers:
            if verbose:
                print(f"  [diagnostico] CLAVE USADA: adata.layers['{k}']")
            return np.asarray(adata.layers[k]), adata, k

    if verbose:
        print("  [diagnostico] NINGUNA clave con shape=n_genes encontrada.")
    return None, adata, None


# --------------------------------------------------------------------------- #
# METRICAS (identicas al resto de scripts de metodos)                         #
# --------------------------------------------------------------------------- #
def evaluar(pred_mat, gt, recalibrar=True, seed=0):
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
def main(indir=INDIR, datasets=DATASETS, csv_out=CSV_OUT):
    resultados = []
    print(f"[datos] {indir}")
    print(f"{'dataset':<20}{'clave':<16}{'RMSE':>9}{'MAE':>9}{'bias':>9}{'pearson':>10}"
          f"{'spearman':>11}{'r2':>9}   |{'r2_recal':>10}")
    print("-" * 114)
    for nombre in datasets:
        # ruta_dataset/ruta_ground_truth resuelven el layout: datasets/ es plano
        # y datasets_gen3/ tiene una subcarpeta por base (ver rutas.py)
        h5 = ruta_dataset(nombre, indir)
        gtf = ruta_ground_truth(nombre, indir)
        if not os.path.exists(h5):
            print(f"  [saltado] no existe {h5}")
            continue
        a = ad.read_h5ad(h5)
        gt = dict(np.load(gtf))
        if len(gt["true_value"]) == 0:
            # el .npz existe pero vacio -> el dataset se genero con mask_frac=0
            # y no hay held-out contra el que medir nada
            print(f"  [saltado] {nombre}: sin held-out (mask_frac=0 en la "
                  f"generacion); regenera con mask_frac>0")
            continue

        pred_mat, _, clave = impute_graphst(a)
        if pred_mat is None:
            print(f"{nombre:<20}  [sin reconstruccion en espacio de genes -- ver diagnostico arriba]")
            continue
        m = evaluar(pred_mat, gt)
        m2 = dict(dataset=nombre, metodo=f"GraphST[{clave}]", **m)
        resultados.append(m2)
        print(f"{nombre:<20}{clave:<16}{m['rmse']:>9.3f}{m['mae']:>9.3f}{m['bias']:>9.3f}"
              f"{m['pearson']:>10.4f}{m['spearman']:>11.4f}{m['r2']:>9.4f}   |"
              f"{m['r2_recal']:>10.4f}")

    if resultados:
        guardar_csv(resultados, path=csv_out)
        print(f"\n[save] {csv_out}  (+{len(resultados)} filas)")
    else:
        print("\nNo se guardo nada: GraphST no produjo reconstruccion evaluable "
              "en ningun dataset (ver diagnosticos arriba).")


def main_barrido():
    """
    Barrido combinado: n_neighbors x dim_output x alpha, solo poisson/nb_s00
    (unicos con senal real detectable). Usa r2_recal para aislar "captura de
    senal" de "escala mal calibrada".
    """
    BARRIDO_N_NEIGHBORS = [3, 5, 10, 20]
    BARRIDO_DIM_OUTPUT = [16, 32, 64]
    BARRIDO_ALPHA = [10, 50, 100]
    DATASETS_CON_SENAL = ["poisson", "nb_s00"]
    BETA_FIJO = 1

    print(f"{'dataset':<10}{'n_neigh':>9}{'dim_out':>9}{'alpha':>7}{'clave':<12}"
          f"{'pearson':>10}{'r2_recal':>10}{'bias_recal':>12}")
    print("-" * 90)
    filas = []
    for nombre in DATASETS_CON_SENAL:
        h5 = os.path.join(INDIR, f"dataset_{nombre}.h5ad")
        gtf = os.path.join(INDIR, f"ground_truth_{nombre}.npz")
        a_base = ad.read_h5ad(h5)
        gt = dict(np.load(gtf))

        for nn in BARRIDO_N_NEIGHBORS:
            for do in BARRIDO_DIM_OUTPUT:
                for al in BARRIDO_ALPHA:
                    pred_mat, _, clave = impute_graphst(
                        a_base, n_neighbors=nn, dim_output=do,
                        alpha=al, beta=BETA_FIJO, verbose=False,
                    )
                    if pred_mat is None:
                        print(f"{nombre:<10}{nn:>9}{do:>9}{al:>7}  [sin reconstruccion]")
                        continue
                    m = evaluar(pred_mat, gt)
                    print(f"{nombre:<10}{nn:>9}{do:>9}{al:>7}{clave:<12}"
                          f"{m['pearson']:>10.4f}{m['r2_recal']:>10.4f}{m['bias_recal']:>12.3f}")
                    filas.append(dict(
                        dataset=nombre, n_neighbors=nn, dim_output=do,
                        alpha=al, beta=BETA_FIJO, clave=clave, **m
                    ))

    if filas:
        path = os.path.join(INDIR, "barrido_graphst_combinado.csv")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
            w.writeheader(); w.writerows(filas)
        print(f"\n[save] {path}")
    else:
        print("\nNingun punto del barrido produjo reconstruccion evaluable.")


if __name__ == "__main__":
    if "--barrido" in sys.argv:
        main_barrido()
    elif "--gen3" in sys.argv:
        main(indir=GEN3, datasets=DATASETS_GEN3, csv_out=CSV_OUT_GEN3)
    else:
        main()
