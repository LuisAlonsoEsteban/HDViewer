"""
remask.py
=========
Cambia el mask_frac de los datasets de UNA base SIN regenerarlos.

Por que re-enmascarar en vez de regenerar: el held-out es una operacion
POSTERIOR a la generacion -- mask_ground_truth() se aplica sobre los counts ya
muestreados y solo decide que entradas se ponen a 0. Los counts subyacentes
(layers["ground_truth_full"]), el campo lambda, los dominios y las coordenadas
no dependen de mask_frac. Regenerar los tocaria todos y, ademas, el fuente de
gen_datasets3.py se perdio y lo unico que hay es una decompilacion (ver
gen_datasets3_recuperado.py): no conviene apoyarse en ella para reconstruir el
ground truth del TFM cuando la alternativa es exacta.

Se reutiliza la mask_ground_truth() REAL de gen_datasets2, la misma que uso el
generador original, asi que el procedimiento de enmascarado es identico.

POR QUE 0.30 EN nb_s00: con mask_frac=0.10 el dataset esta saturado y no
discrimina. Medido sobre 5 semillas con la rejilla de resolucion actual
(techo = ARI de ground_truth_full, que no depende de mask_frac):

    mask_frac   suelo nb_s00        margen (techo-suelo)
      0.10      0.995 +- 0.007          +0.004
      0.20      0.923 +- 0.041          +0.076
      0.30      0.726 +- 0.077          +0.273
      0.40      0.387 +- 0.229          +0.612   <- sd se dispara
      0.50      0.052 +- 0.104          +0.947   <- suelo degenerado

0.30 es el mayor margen con el suelo aun estable. nb_s10 se queda en 0.10: alli
0.10 ya da margen +0.227 y subirlo degenera el suelo a 0.000 en 0.40.

El mask_frac pasa a ser POR BASE, no global. Es defendible porque las dos bases
parten de sparsity nativa distinta (nb_s00 ~0%, nb_s10 ~10%), asi que enmascarar
la misma fraccion no supone la misma perdida efectiva de informacion.

Uso:
    python remask.py --base nb_s00 --mask-frac 0.30
    python remask.py --base nb_s00 --mask-frac 0.30 --dry-run
"""

import os
import sys
import glob
import zlib
import argparse
import numpy as np
import anndata as ad

from rutas import GEN3, dir_base, ruta_dataset, ruta_ground_truth
from gen_datasets2 import mask_ground_truth


def semilla_de(nombre):
    """Semilla determinista y estable entre ejecuciones derivada del nombre del
    dataset. No se usa hash() de Python: su valor para str cambia entre procesos
    (PYTHONHASHSEED), asi que re-enmascarar dos veces daria mascaras distintas."""
    return zlib.crc32(nombre.encode("utf-8"))


def remask_uno(nombre, mask_frac, raiz, dry_run=False):
    ruta_h5 = ruta_dataset(nombre, raiz)
    a = ad.read_h5ad(ruta_h5)

    full = np.asarray(a.layers["ground_truth_full"], dtype=np.float64)
    rng = np.random.default_rng(semilla_de(nombre))
    obs, mask, rows, cols, valores = mask_ground_truth(full, mask_frac, rng)

    mf_viejo = a.uns["generation_params"].get("mask_frac")
    sp_antes = 100.0 * float((np.asarray(a.X) == 0).mean())
    sp_ahora = 100.0 * float((obs == 0).mean())

    if dry_run:
        print(f"  [dry-run] {nombre:26} mask_frac {mf_viejo} -> {mask_frac}  "
              f"held-out {int(a.layers['mask'].sum())} -> {int(mask.sum())}  "
              f"sparsity {sp_antes:.2f}% -> {sp_ahora:.2f}%")
        return

    a.X = obs.astype(np.float32)
    a.layers["mask"] = mask.astype(np.uint8)
    a.uns["generation_params"]["mask_frac"] = float(mask_frac)
    # rastro de la operacion: sin esto, un .h5ad re-enmascarado es
    # indistinguible de uno generado con ese mask_frac desde el principio
    a.uns["generation_params"]["mask_frac_original"] = float(mf_viejo)
    a.uns["generation_params"]["remask_semilla"] = int(semilla_de(nombre))

    a.write_h5ad(ruta_h5)
    np.savez(ruta_ground_truth(nombre, raiz),
             spot_idx=rows.astype(np.int64),
             gene_idx=cols.astype(np.int64),
             true_value=valores.astype(np.float64))

    print(f"  {nombre:26} mask_frac {mf_viejo} -> {mask_frac}  "
          f"held-out {int(mask.sum())}  sparsity {sp_antes:.2f}% -> {sp_ahora:.2f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="base a re-enmascarar (p.ej. nb_s00)")
    ap.add_argument("--mask-frac", type=float, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    carpeta = dir_base(args.base, GEN3)
    nombres = sorted(os.path.basename(f)[len("dataset_"):-len(".h5ad")]
                     for f in glob.glob(os.path.join(carpeta, "dataset_*.h5ad")))
    if not nombres:
        raise SystemExit(f"No hay datasets de {args.base} en {carpeta}")

    print(f"Re-enmascarando {len(nombres)} datasets de {args.base} "
          f"a mask_frac={args.mask_frac}\n")
    for n in nombres:
        remask_uno(n, args.mask_frac, GEN3, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
