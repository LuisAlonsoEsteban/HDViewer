"""
combinar_benchmarks.py
======================
Une los CSV por metodo que produce evaluar_clustering_metodos.py
(benchmark_*_dominios_<metodo>.csv) en las dos tablas finales.

Las filas de referencia sin_heldout / sin_imputar las calcula cada job por
separado, asi que aparecen repetidas: se deduplican por (dataset, metodo).

Uso:
    python combinar_benchmarks.py
"""

import os
import glob
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


def combinar(patron, salida):
    # Se excluyen las salidas DERIVADAS, no solo `salida`: el patron
    # "benchmark_clustering_dominios_*.csv" tambien casa con
    # "..._agregado.csv", que lo produce agregar_semillas.py a partir de esta
    # misma tabla. Colarlo aqui metia filas con dataset=NaN y columnas
    # ARI_mean/ARI_std que no existen en los CSV por metodo.
    DERIVADOS = ("_agregado.csv",)
    ficheros = sorted(f for f in glob.glob(os.path.join(OUTDIR, patron))
                      if not f.endswith(salida)
                      and not f.endswith(DERIVADOS))
    if not ficheros:
        print(f"[aviso] sin ficheros para {patron}")
        return None
    df = pd.concat([pd.read_csv(f) for f in ficheros], ignore_index=True)
    df = df.drop_duplicates(subset=["dataset", "metodo"], keep="first")
    df["_o"] = df["metodo"].apply(clave_orden)
    df = df.sort_values(["dataset", "_o"]).drop(columns="_o")
    ruta = os.path.join(OUTDIR, salida)
    df.to_csv(ruta, index=False)
    print(f"\n{ruta}   ({len(ficheros)} ficheros combinados)")
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    combinar("benchmark_clustering_dominios_*.csv",
             "benchmark_clustering_dominios.csv")
    combinar("benchmark_imputacion_dominios_*.csv",
             "benchmark_imputacion_dominios.csv")
