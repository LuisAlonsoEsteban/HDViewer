"""
evaluar_clustering_metodos.py
==============================
Generaliza evaluar_clustering_knn.py a los 4 metodos de imputacion
VALIDOS (kNN, SpaMean, STAGATE, MAGIC). GraphST queda excluido: ya
documentado que su obsm['emb'] es un embedding de autoencoder, no una
reconstruccion en espacio de genes (bias~-46 en las 4 bases, ~valor de
la media, y r2 fuertemente negativo -- no es comparable pese a que
dim_output=100 coincida por casualidad con n_genes).

Para cada dataset x metodo:
    - imputa la matriz enmascarada
    - mide RMSE/R2 sobre las posiciones held-out (igual que benchmark_metodos.csv)
    - corre Leiden (cluster_domains.py) sobre la matriz imputada
    - mide ARI/NMI/AMI/HOMO contra el dominio real
    - reporta frac_recuperada = ARI_metodo / ARI_techo(sin_heldout)

Reutiliza los adaptadores YA VALIDADOS en time_methods.py.

Mide ademas, POR METODO, el coste de la llamada de imputacion en si:
t_wall_s / t_cpu_s y el pico de RSS (rss_pico_mb, rss_delta_mb). Esas
columnas van en ambos CSV.

Salidas en datasets_gen3/, con sufijo por metodo para que los jobs en
paralelo no se pisen el fichero:
    benchmark_clustering_dominios_<metodo>.csv  (ARI/NMI/AMI/HOMO + recursos)
    benchmark_imputacion_dominios_<metodo>.csv  (RMSE/R2 held-out + recursos)

Uso:
    python evaluar_clustering_metodos.py --metodo knn
    python evaluar_clustering_metodos.py --metodo all      # los cuatro
"""

import os
import sys
import time
import argparse
import threading
import warnings
import numpy as np
import pandas as pd
import psutil
from scipy import stats

from cluster_domains import (OUTDIR, CFG, cargar, construir_grafo,
                            elegir_resolucion, evaluar, datasets_disponibles)
from rutas import METHODS, ruta_ground_truth
sys.path.insert(0, METHODS)

# adaptadores YA VALIDADOS, reutilizados tal cual de time_methods.py
from time_methods import _run_knn, _run_spamean, _run_stagate, _run_magic

# Por defecto TODAS las variantes "libre" que haya en datasets_gen3/ (una por
# cada base de gen_datasets2: gp_crudo, gaussian, poisson, nb_s00, nb_s10,
# nb_s50), descubiertas con datasets_disponibles() -- antes eran solo las dos
# NB porque eran las unicas que gen_datasets3 generaba. Las "atada" quedan
# fuera del default: "libre" es el modo que build_domain_offsets documenta como
# recomendado (aisla la senal de dominio de la de modulo); se anaden con
# --datasets. Los jobs de SLURM pasan su propia lista (las replicas _seedN), asi
# que este default solo afecta a las ejecuciones locales.
DATASETS = datasets_disponibles(patron="_dom_libre")

# --------------------------------------------------------------------------- #
# PREPROCESADO DEL CLUSTERING: DECISION EXPLICITA, NO INFERIDA                 #
# --------------------------------------------------------------------------- #
# Antes se dejaba en modo AUTO (construir_grafo -> es_matriz_de_conteo), y eso
# hacia que la rama de preprocesado dependiera de la SALIDA de cada metodo: un
# metodo que devolviera un solo valor negativo (un decoder lineal, una
# factorizacion, cualquier cosa sin ReLU final) caia en la otra rama sin avisar
# y su ARI dejaba de ser comparable con el de los demas.
#
# Se fija a False para TODOS. Motivo medido (nb_s00/nb_s10 _dom_libre_seed42,
# forzando las dos ramas sobre las mismas matrices):
#
#     matriz                normalizar=True     normalizar=False
#     nb_s00 techo               0.961 (4c)         0.958 (4c)
#     nb_s00 suelo               0.000 (1c)         1.000 (3c)
#     nb_s10 techo               0.559 (2c)         0.842 (3c)
#     nb_s10 suelo               0.000 (1c)         0.617 (3c)
#
# Con normalizar=True el suelo colapsa a 1 cluster y el "margen de imputacion"
# que abre el held-out es un artefacto del preprocesado, no una propiedad del
# dato: en nb_s00 el margen real es CERO (suelo 1.000 >= techo 0.958). El
# mecanismo no es el tamano de libreria (medido: identico entre dominios,
# 4904/4915/4858, y la fraccion enmascarada es un 10% uniforme), sino la
# interaccion de log1p+scale con el 10% de ceros exactos de la mascara.
#
# Los conteos de gen_datasets3 ya salen con tamano de libreria controlado, asi
# que normalize_total no esta corrigiendo nada real aqui.
NORMALIZAR = False

# Un fichero POR METODO: cuando los 4 jobs corren en paralelo, un unico CSV
# compartido se pisa entre ellos (solo sobrevive el ultimo en terminar).
# combinar_csv.py los une despues en la tabla final.
def rutas_csv(sufijo):
    return (os.path.join(OUTDIR, f"benchmark_clustering_dominios_{sufijo}.csv"),
            os.path.join(OUTDIR, f"benchmark_imputacion_dominios_{sufijo}.csv"))

# GraphST excluido: obsm['emb'] es un embedding de autoencoder, no
# reconstruccion en espacio de genes (ver benchmark_metodos.csv, bias~-46
# en las 4 bases -- documentado, no un descarte improvisado aqui)
METODOS = {
    "knn":     ("kNN_espacial(k=8)", _run_knn),
    "spamean": ("SpaMean",           _run_spamean),
    "stagate": ("STAGATE",           _run_stagate),
    "magic":   ("MAGIC",             _run_magic),
}


def medir(run_fn, a, intervalo=0.05):
    """Ejecuta run_fn(a) midiendo SOLO esa llamada: tiempo de pared, tiempo de
    CPU (suma de todos los hilos: MAGIC y STAGATE paralelizan) y pico de RSS.

    El pico se muestrea con un hilo aparte en vez de leer ru_maxrss porque
    ru_maxrss es el maximo del proceso ENTERO desde que arranco -- una vez que
    TensorFlow reserva sus buffers ya no baja, asi que daria el mismo valor a
    los cuatro metodos. Lo que interesa es rss_delta_mb: cuanto sube el RSS
    por encima de la linea base inmediatamente anterior a la llamada, que si
    es atribuible al metodo."""
    proc = psutil.Process()
    base = proc.memory_info().rss
    pico = [base]
    parar = threading.Event()

    def muestrear():
        while not parar.is_set():
            try:
                pico[0] = max(pico[0], proc.memory_info().rss)
            except psutil.Error:
                pass
            parar.wait(intervalo)

    hilo = threading.Thread(target=muestrear, daemon=True)
    hilo.start()
    t0, c0 = time.perf_counter(), time.process_time()
    try:
        pred = run_fn(a)
    finally:
        parar.set()
        hilo.join()
    t_wall = time.perf_counter() - t0
    t_cpu = time.process_time() - c0

    MB = 1024 ** 2
    return pred, dict(t_wall_s=round(t_wall, 3), t_cpu_s=round(t_cpu, 3),
                      rss_base_mb=round(base / MB, 1),
                      rss_pico_mb=round(pico[0] / MB, 1),
                      rss_delta_mb=round((pico[0] - base) / MB, 1))


def calc_metricas_imputacion(pred_mat, gt):
    """Identico patron a evaluar() en impute_knn_spatial.py: RMSE/MAE/bias/
    pearson/spearman/r2 + recalibracion afin, sobre las posiciones held-out."""
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

    # Diagnostico de la SALIDA CRUDA del metodo, no de las posiciones held-out.
    # Es la comprobacion que decide si el metodo es compatible con un
    # preprocesado de conteo: si pred_neg_pct > 0, normalize_total/log1p sobre
    # esa matriz daria NaN. Va al CSV para que quede registrado por metodo en
    # vez de tener que reproducirlo a mano.
    todo = np.asarray(pred_mat, dtype=float)
    return dict(n_heldout=len(true), rmse=rmse, mae=mae, bias=bias,
               pearson=pearson, spearman=spearman, r2=r2, r2_recal=float(r2_recal),
               pred_min=float(np.nanmin(todo)), pred_max=float(np.nanmax(todo)),
               pred_neg_pct=round(100.0 * float((todo < 0).mean()), 4),
               pred_nan=int(np.isnan(todo).sum()))


def correr_clustering(a, dom_real, n_target, cfg, etiqueta, nombre_dataset,
                      normalizar=NORMALIZAR, modo="referencia"):
    """normalizar se pasa SIEMPRE explicito (ver NORMALIZAR arriba): todas las
    filas de la tabla -- techo, suelo y los N metodos -- tienen que pasar por el
    mismo preprocesado o el ARI no es comparable entre ellas."""
    b = construir_grafo(a, cfg, normalizar=normalizar)
    res, leiden, nclust, barrido = elegir_resolucion(b, n_target, cfg)
    m = evaluar(dom_real, leiden)
    # El barrido COMPLETO va al CSV, no solo la resolucion ganadora. Sin el, un
    # n_clusters != n_dominios es ambiguo: no se distingue "ninguna resolucion
    # del grid daba n_dominios" (grid insuficiente -> se arregla ampliando el
    # grid) de "alguna lo daba y el criterio eligio otra" (criterio de seleccion
    # -> se arregla en elegir_resolucion). Son dos arreglos distintos, y el
    # primero era justo lo que hundia el techo de nb_s00 por debajo del suelo.
    m.update(dataset=nombre_dataset, metodo=etiqueta, modo=modo,
             n_dominios=n_target,
             normalizado=bool(b.uns["normalizado"]),
             barrido=" ".join(f"r{r}->{n}c" for r, n in barrido),
             sparsity_pct=round(100 * float((np.asarray(a.X) == 0).mean()), 2),
             resolution=res)
    print(f"[{nombre_dataset:>18} | {etiqueta:>18}]  norm={str(m['normalizado']):>5}  "
          f"sparsity={m['sparsity_pct']:>5.2f}%  "
          f"res={res}  n_clust={m['n_clusters']:>2}  ARI={m['ARI']:.3f}  HOMO={m['HOMO']:.3f}")
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metodo", default="all",
                    choices=list(METODOS.keys()) + ["all"],
                    help="metodo a evaluar; 'all' corre los cuatro en el mismo job")
    ap.add_argument("--datasets", nargs="+", default=DATASETS,
                    help="datasets de datasets_gen3/ a evaluar")
    ap.add_argument("--normalizar", choices=["si", "no"],
                    default="si" if NORMALIZAR else "no",
                    help="normalize_total+log1p antes del PCA. Se aplica por "
                         "IGUAL a techo, suelo y metodos. Por defecto 'no' "
                         "(ver NORMALIZAR en la cabecera del fichero)")
    args = ap.parse_args()
    normalizar = (args.normalizar == "si")

    seleccion = (METODOS if args.metodo == "all"
                 else {args.metodo: METODOS[args.metodo]})
    csv_cluster, csv_imput = rutas_csv(args.metodo)

    # NO se silencia todo: el aviso de construir_grafo() cuando se fuerza
    # normalizar=True sobre una matriz con negativos es un RuntimeWarning, y con
    # filterwarnings("ignore") a secas quedaria enterrado justo el mensaje que
    # avisa de que el metodo cayo fuera del pipeline comun.
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    filas_cluster, filas_imput = [], []

    print(f"Metodos a evaluar: {', '.join(n for n, _ in seleccion.values())}")
    print(f"Datasets: {', '.join(args.datasets)}")
    print(f"Preprocesado del clustering: normalize_total+log1p = {normalizar} "
          f"(mismo para techo, suelo y todos los metodos)\n")

    for nombre in args.datasets:
        a = cargar(nombre)
        dom_real = a.obs["dominio"].astype(int).values
        n_target = int(len(np.unique(dom_real)))
        gt = np.load(ruta_ground_truth(nombre, OUTDIR))

        # --- techo: sin enmascarar (referencia ideal, igual para todos los metodos) ---
        a_full = a.copy()
        a_full.X = np.asarray(a.layers["ground_truth_full"], dtype=np.float32)
        m_techo = correr_clustering(a_full, dom_real, n_target, CFG,
                                    "sin_heldout", nombre, normalizar)
        filas_cluster.append(m_techo)

        # --- suelo: observado con mascara, sin imputar ---
        m_suelo = correr_clustering(a, dom_real, n_target, CFG,
                                    "sin_imputar", nombre, normalizar)
        filas_cluster.append(m_suelo)

        for metodo_nombre, run_fn in seleccion.values():
            pred, recursos = medir(run_fn, a)
            print(f"[{nombre:>18} | {metodo_nombre:>18}]  "
                  f"t_wall={recursos['t_wall_s']:>8.3f}s  "
                  f"t_cpu={recursos['t_cpu_s']:>8.3f}s  "
                  f"RSS pico={recursos['rss_pico_mb']:>7.1f}MB  "
                  f"(+{recursos['rss_delta_mb']:.1f}MB sobre base)")

            m_imput = calc_metricas_imputacion(pred, gt)
            m_imput.update(dataset=nombre, metodo=metodo_nombre, **recursos)
            filas_imput.append(m_imput)

            # --- DENOISE vs IMPUTA ------------------------------------------
            # Antes se evaluaba UNA sola matriz, `pred` entera, y eso mezclaba
            # dos tareas distintas sin decirlo. Sustituir la matriz completa
            # deja que el metodo reescriba tambien las entradas OBSERVADAS: eso
            # es denoising. Restaurar los observados y dejar solo los held-out
            # imputados es imputacion. Los cuatro metodos no se comportaban
            # igual respecto a eso (SpaMean conservaba ceros, kNN/MAGIC/STAGATE
            # reescribian las 40.000 entradas), asi que la comparacion no era
            # homogenea y el eje en el que ganaba cada uno era invisible.
            #
            # Se reportan los dos, con la misma imputacion detras: la diferencia
            # entre las dos filas es exactamente cuanto de su ARI viene de
            # suavizar el dato observado y cuanto de rellenar el hueco.
            X_obs = np.asarray(a.X, dtype=np.float32)
            mascara = np.zeros(X_obs.shape, dtype=bool)
            mascara[gt["spot_idx"], gt["gene_idx"]] = True

            for modo, X_eval in (
                ("denoise", pred.astype(np.float32)),
                ("imputa",  np.where(mascara, pred, X_obs).astype(np.float32)),
            ):
                a_eval = a.copy()
                a_eval.X = X_eval
                m_clust = correr_clustering(a_eval, dom_real, n_target, CFG,
                                            f"{metodo_nombre}[{modo}]", nombre,
                                            normalizar, modo=modo)
                # ARI ABSOLUTO como metrica principal, con techo y suelo al lado.
                # frac_recuperada ahora normaliza al intervalo [suelo, techo] en
                # vez de dividir por el techo: con el suelo honesto, dividir por
                # el techo mentia de forma visible -- kNN en nb_s10 daba
                # frac_recuperada=0.596 con un ARI (0.432) POR DEBAJO del suelo
                # (0.550), o sea "recupera el 60%" de algo que en realidad
                # destruye. Con el anclaje a los dos extremos, un metodo que
                # empeora el suelo sale negativo, que es lo que se quiere leer.
                #
                # OJO: el techo NO es un maximo estricto -- un metodo que suaviza
                # puede superarlo y dar >1. Metrica secundaria; el ARI crudo con
                # sus dos anclajes al lado es lo que se reporta.
                m_clust["ARI_techo"] = m_techo["ARI"]
                m_clust["ARI_suelo"] = m_suelo["ARI"]
                m_clust["frac_recuperada"] = (
                    (m_clust["ARI"] - m_suelo["ARI"]) /
                    (m_techo["ARI"] - m_suelo["ARI"])
                    if m_techo["ARI"] > m_suelo["ARI"] else float("nan"))
                m_clust.update(recursos)
                filas_cluster.append(m_clust)

    cols_rec = ["t_wall_s", "t_cpu_s", "rss_base_mb", "rss_pico_mb", "rss_delta_mb"]
    cols_c = ["dataset", "metodo", "modo", "normalizado", "sparsity_pct", "n_dominios",
              "n_clusters", "ARI", "ARI_techo", "ARI_suelo", "NMI", "AMI", "HOMO",
              "resolution", "barrido", "frac_recuperada"] + cols_rec
    df_c = pd.DataFrame(filas_cluster)
    for c in cols_c:
        if c not in df_c.columns:
            df_c[c] = np.nan
    df_c = df_c[cols_c].round(4)
    df_c.to_csv(csv_cluster, index=False)

    cols_i = ["dataset", "metodo", "n_heldout", "rmse", "mae", "bias",
              "pearson", "spearman", "r2", "r2_recal",
              "pred_min", "pred_max", "pred_neg_pct", "pred_nan"] + cols_rec
    df_i = pd.DataFrame(filas_imput)[cols_i].round(4)
    df_i.to_csv(csv_imput, index=False)

    print(f"\nTabla clustering -> {csv_cluster}")
    print(df_c.to_string(index=False))
    print(f"\nTabla imputacion (held-out) -> {csv_imput}")
    print(df_i.to_string(index=False))


if __name__ == "__main__":
    main()