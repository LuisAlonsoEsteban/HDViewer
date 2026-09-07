"""
agregar_semillas.py
===================
Agrega las replicas por semilla del benchmark: media +- sd de cada metrica
sobre las N semillas, por (dataset_base, metodo).

Una replica = una realizacion completa del proceso generativo con otra semilla
(campo GP, particion de dominios, offsets y counts). Con una sola semilla no se
puede saber si la diferencia entre dos metodos es real o ruido de esa tirada
concreta: el techo de nb_s10_dom_libre se movio 0.559 -> 0.736 solo por cambiar
la semilla, un margen del orden de las diferencias entre metodos.

Lee benchmark_clustering_dominios.csv / benchmark_imputacion_dominios.csv (los
que produce combinar_benchmarks.py), donde la columna 'dataset' trae la semilla
en el nombre (nb_s10_dom_libre_seed7).

Uso:
    python combinar_benchmarks.py && python agregar_semillas.py
"""

import os
import re
import numpy as np
import pandas as pd

from cluster_domains import OUTDIR

ORDEN = ["sin_heldout", "sin_imputar", "kNN_espacial(k=8)", "SpaMean",
         "STAGATE", "MAGIC"]


def clave_orden(metodo):
    """Posicion en ORDEN ignorando el sufijo de modo: evaluar_clustering_metodos
    emite "MAGIC[denoise]" / "MAGIC[imputa]", que con una busqueda literal caian
    todos al final en bloque en vez de junto a su metodo."""
    base = str(metodo).split("[", 1)[0]
    return ORDEN.index(base) if base in ORDEN else len(ORDEN)


def separar_semilla(df):
    ex = df["dataset"].str.extract(r"^(?P<base>.+?)_seed(?P<seed>\d+)$")
    df = df.copy()
    df["base"] = ex["base"].fillna(df["dataset"])
    df["seed"] = ex["seed"]
    return df


def agregar(nombre_csv, metricas, salida):
    ruta = os.path.join(OUTDIR, nombre_csv)
    if not os.path.exists(ruta):
        print(f"[aviso] falta {ruta}")
        return
    df = separar_semilla(pd.read_csv(ruta))
    n_seeds = df["seed"].nunique(dropna=True)
    metricas = [m for m in metricas if m in df.columns]

    g = df.groupby(["base", "metodo"])[metricas]
    out = g.agg(["mean", "std", "count"])
    out.columns = [f"{m}_{s}" for m, s in out.columns]
    out = out.reset_index()
    out["_o"] = out["metodo"].apply(clave_orden)
    out = out.sort_values(["base", "_o"]).drop(columns="_o")
    out.to_csv(os.path.join(OUTDIR, salida), index=False)

    print(f"\n=== {nombre_csv}  ({n_seeds} semillas) -> {salida}")
    for base, sub in out.groupby("base"):
        print(f"\n  {base}")
        for _, r in sub.iterrows():
            partes = []
            for m in metricas:
                mu, sd, n = r[f"{m}_mean"], r[f"{m}_std"], r[f"{m}_count"]
                if np.isnan(mu):
                    continue
                sd = 0.0 if np.isnan(sd) else sd
                partes.append(f"{m}={mu:6.3f}+-{sd:5.3f}")
            print(f"    {r['metodo']:>18}  n={int(r[metricas[0]+'_count']):>2}  "
                  + "  ".join(partes))
    return out


if __name__ == "__main__":
    agregar("benchmark_clustering_dominios.csv",
            ["ARI", "NMI", "HOMO", "frac_recuperada", "t_wall_s", "rss_delta_mb"],
            "benchmark_clustering_dominios_agregado.csv")
    agregar("benchmark_imputacion_dominios.csv",
            ["rmse", "pearson", "r2", "r2_recal", "t_wall_s"],
            "benchmark_imputacion_dominios_agregado.csv")
