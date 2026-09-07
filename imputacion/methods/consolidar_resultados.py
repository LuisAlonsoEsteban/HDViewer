"""
consolidar_resultados.py
==========================
Lee benchmark_metodos.csv (kNN, SpaMean, STAGATE, GraphST, MAGIC) y genera
la tabla final comparativa usando r2_recal (recalibracion afin honesta,
holdout) para los 5 metodos con el mismo criterio.

GraphST se incluye pero se anota como "no aplicable" en la conclusion
(Pearson≈0.01-0.03 en barrido de 72 configs, indistinguible de ruido).

Uso:
    python consolidar_resultados.py
"""

import pandas as pd
import os

# rutas.py vive en imputacion/, el directorio padre de methods/: se anade al
# path para importarlo sin depender del CWD ni de rutas escritas a mano
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rutas import DATASETS as INDIR
BENCHMARK_CSV = os.path.join(INDIR, "benchmark_metodos.csv")
DATASETS = ["poisson", "nb_s00", "nb_s10", "nb_s50"]

METODOS = {
    "kNN espacial":   "kNN_espacial(k=8)",
    "SpaMean-Impute": "SpaMean-Impute(k=8,delta=0.4)",
    "STAGATE_pyG":    "STAGATE_pyG_crudo",
    "GraphST":        "GraphST",
    "MAGIC":          "MAGIC",
}


def main():
    df = pd.read_csv(BENCHMARK_CSV)
    df = df.drop_duplicates(subset=["dataset", "metodo"], keep="last")

    filas = []
    for nombre in DATASETS:
        fila = {"dataset": nombre}
        for metodo_final, metodo_csv in METODOS.items():
            sub = df[(df["dataset"] == nombre) & (df["metodo"] == metodo_csv)]
            fila[metodo_final] = round(sub["r2_recal"].values[0], 4) if len(sub) else None
        filas.append(fila)

    tabla = pd.DataFrame(filas).set_index("dataset")

    print("\nTABLA FINAL -- R2 recalibrado (holdout, todos los metodos)\n")
    print(tabla.to_string())

    out_csv = os.path.join(INDIR, "tabla_final_comparativa.csv")
    tabla.to_csv(out_csv)
    print(f"\n[save] {out_csv}")

    out_md = os.path.join(INDIR, "tabla_final_comparativa.md")
    with open(out_md, "w") as f:
        f.write(tabla.to_markdown())
        f.write(
            "\n\n**Nota sobre GraphST**: la reconstruccion interna "
            "(obsm['emb'], en realidad el tensor `h` del decoder de "
            "Encoder.forward()) fue evaluada tras un barrido combinado de "
            "72 configuraciones (n_neighbors x dim_output x alpha, este "
            "ultimo confirmado sin efecto -- ver bug pendiente). El maximo "
            "Pearson alcanzado en cualquier configuracion fue ~0.03, "
            "indistinguible de ruido. Se reporta en la tabla por "
            "completitud, pero se considera no aplicable a imputacion "
            "gen-a-gen en este dataset sintetico.\n"
        )
    print(f"[save] {out_md}")


if __name__ == "__main__":
    main()