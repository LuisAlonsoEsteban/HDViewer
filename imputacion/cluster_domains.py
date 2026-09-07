"""
cluster_domains.py
==================
Experimento B (recuperacion de dominios): sobre TODOS los datasets con dominios
de gen_datasets3 ({gp_crudo, gaussian, poisson, nb_s00, nb_s10, nb_s50} x
{libre, atada}, descubiertos con datasets_disponibles()), corre un pipeline de
clustering Leiden
sobre la EXPRESION y mide cuanto se parecen los clusters resultantes a los
dominios espaciales reales (obs["dominio"], el ground truth categorico).

NO hay imputacion aqui, pero desde que gen_datasets3 genera con mask_frac=0.10
X ya NO es la matriz honesta: lleva un 10% de genes por spot puestos a 0 como
held-out. Por eso cada dataset se evalua sobre DOS matrices:

    observado  = a.X                        -> el dato degradado, la entrada
                                               real del imputador (SUELO)
    sin_heldout = layers["ground_truth_full"] -> los counts sin enmascarar, con
                                               solo su sparsity biologica
                                               (TECHO alcanzable imputando)

La diferencia de ARI entre ambas filas es exactamente el margen que puede
recuperar un metodo de imputacion en ese dataset.

Metricas (todas via sklearn, ground_truth = dominio real, pred = cluster Leiden):
    ARI  Adjusted Rand Index         - acuerdo de particiones, corregido por azar (1=perfecto, 0=azar)
    NMI  Normalized Mutual Info      - informacion compartida normalizada (1=perfecto)
    AMI  Adjusted Mutual Info        - NMI corregido por azar
    HOMO Homogeneidad                - cada cluster contiene un solo dominio? (1=puro)

Ademas registra n_clusters encontrados vs n_dominios reales (para ver si Leiden
sobre/infra-segmenta) y la resolucion usada.

Salidas en datasets_gen3/:
    metricas_clustering.csv                  - una fila por (dataset, matriz)
    fig_clustering_<dataset>_<matriz>.png    - mapa espacial: dominio real vs
                                               clusters Leiden

Uso:
    python cluster_domains.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    adjusted_rand_score, normalized_mutual_info_score,
    adjusted_mutual_info_score, homogeneity_score,
)

from rutas import (GEN3 as OUTDIR, ruta_dataset, datasets_en,
                   dir_base, base_de)


def datasets_disponibles(indir=None, patron="_dom_", con_semillas=False):
    """Los datasets con dominios que hay realmente en la carpeta.

    Se descubren en vez de listarse a mano porque gen_datasets3.py ya no genera
    solo nb_s00 y nb_s10: genera las SEIS bases de gen_datasets2 (gp_crudo,
    gaussian, poisson, nb_s00, nb_s10, nb_s50) x {libre, atada}. Con la lista
    fija, anadir una base obligaba a tocar tambien este fichero y era facil
    dejarse datasets sin ARI.

    con_semillas=False excluye las replicas <nombre>_seed<N> del barrido: son
    otra pregunta (variabilidad entre realizaciones, ver agregar_semillas.py) y
    multiplicarian por N el coste de esta tabla.
    """
    import re
    indir = indir or OUTDIR
    # datasets_en recorre las subcarpetas por base (ver rutas.py): el layout de
    # la carpeta no se replica aqui
    nombres = [n for n in datasets_en(indir) if patron in f"dataset_{n}.h5ad"]
    if not con_semillas:
        nombres = [n for n in nombres if not re.search(r"_seed\d+$", n)]
    # orden estable por verosimilitud (el de gen_datasets2) y luego por modo,
    # para que la tabla se lea de menos a mas ruido/sparsity
    orden_base = ["gp_crudo", "gaussian", "poisson", "nb_s00", "nb_s10", "nb_s50"]
    def clave(n):
        base, _, modo = n.partition("_dom_")
        b = orden_base.index(base) if base in orden_base else len(orden_base)
        return (b, base, modo != "libre", modo)
    return sorted(nombres, key=clave)


DATASETS = datasets_disponibles()

# --- parametros del pipeline de clustering (fijos para todos los datasets, para
#     que la comparacion nb_s00 vs nb_s10 sea limpia: solo cambia el dato) ---
CFG = dict(
    seed          = 0,
    n_pcs         = 30,
    n_neighbors   = 15,
    # normalize_total + log1p (pipeline scanpy estandar). Se aplica SOLO a los
    # datasets de conteo: ver es_matriz_de_conteo(). Con False se desactiva
    # para todos.
    normalizar    = True,
    # Grid de resolucion: se elige la que da n_clusters mas cercano al
    # n_dominios real y se reporta el ARI ahi (desempate -> resolucion mas
    # baja). Es el enfoque estandar en benchmarking de clustering espacial;
    # evita el artefacto de sobre-segmentacion, no busca un ARI "optimo".
    #
    # AMPLIADO POR ABAJO (antes empezaba en 0.3). El grid viejo no contenia
    # NINGUNA resolucion que diera 3 clusters en el techo de nb_s00, y eso
    # hundia el techo por debajo del suelo -- que es lo que hacia que el
    # "margen de imputacion" saliera negativo. Medido sobre el techo de
    # nb_s00_dom_libre, barriendo r de 0.05 a 1.5:
    #
    #   seed42  r0.05-0.2 -> 3c/ARI 1.00 | r0.3 -> 4c/ARI 0.96
    #   seed13  r0.05-0.2 -> 3c/ARI 1.00 | r0.3 -> 4c/ARI 0.70
    #   seed7   r0.05-0.3 -> 3c/ARI 0.99
    #
    # El tramo alto se conserva entero porque nb_s10 lo necesita en el otro
    # sentido: ahi todo r<=0.3 colapsa a 1 cluster y el optimo esta en 0.5-0.8.
    # Un solo grid tiene que cubrir los dos regimenes.
    #
    # Coste: 9 corridas de Leiden por matriz en vez de 5. Leiden sobre 400
    # spots son ~0.1 s; despreciable frente a los ~150 s de STAGATE.
    res_grid      = (0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5),
)


# matrices sobre las que se corre el MISMO pipeline: etiqueta -> capa del h5ad
# (None = a.X). "observado" es el dato con held-out (suelo, entrada del
# imputador); "sin_heldout" son los counts sin enmascarar (techo).
MATRICES = [("observado", None), ("sin_heldout", "ground_truth_full")]


def cargar(nombre, capa=None):
    """Carga el dataset. Si capa is not None, sustituye X por esa layer -- asi
    todo el pipeline aguas abajo es literalmente el mismo codigo y lo unico que
    cambia entre filas es la matriz de expresion."""
    a = ad.read_h5ad(ruta_dataset(nombre, OUTDIR))
    # ground truth categorico de dominio
    a.obs["dominio"] = a.obs["dominio"].astype("category")
    if capa is not None:
        a.X = np.asarray(a.layers[capa], dtype=np.float32)
    return a


def es_matriz_de_conteo(a):
    """Decide si a esta matriz le toca la normalizacion estandar de scanpy
    (normalize_total + log1p) o no.

    Esa transformacion PRESUPONE conteos no negativos. Aplicarla a las
    verosimilitudes continuas de gen_datasets3 es un error de tipo: con
    'gaussian' hay valores negativos (~0.1% con mu=50), log1p los convierte en
    NaN y el PCA revienta; con 'gp_crudo' no revienta pero se estaria
    log-transformando un campo que no es un conteo.

    Dos condiciones, no una:
      - la VEROSIMILITUD del dataset es de conteo (poisson / nbinom), leida de
        generation_params (gen_datasets3 la guarda explicita);
      - y la matriz que se va a usar no tiene negativos -- importa porque este
        mismo pipeline se corre sobre matrices IMPUTADAS
        (evaluar_clustering_metodos.py), y un metodo puede devolver negativos
        sobre un dataset de conteo.
    """
    p = dict(a.uns.get("generation_params", {}))
    veros = str(p.get("verosimilitud", p.get("modelo", "")))
    es_conteo = veros in ("poisson", "nbinom") or veros.startswith("nb")
    X = np.asarray(a.X, dtype=float)
    return bool(es_conteo and np.nanmin(X) >= 0.0)


def construir_grafo(a, cfg, normalizar=None):
    """Pipeline scanpy hasta el grafo de vecinos (la parte cara), UNA vez por
    dataset. Leiden luego se corre barato en cada resolucion del grid.

    normalizar=None  -> se decide con es_matriz_de_conteo() (modo AUTO, el que
                        usa cluster_domains.main()).
    normalizar=bool  -> se respeta tal cual. es_matriz_de_conteo() sigue
                        corriendo, pero solo para AVISAR, no para cambiar de
                        rama.

    Por que existe la distincion: en modo AUTO la decision depende del CONTENIDO
    de la matriz, y este mismo pipeline se corre sobre matrices IMPUTADAS
    (evaluar_clustering_metodos.py). Un metodo que devolviera un solo valor
    negativo caia en la otra rama de preprocesado sin avisar, y la rama no es un
    detalle: medido sobre nb_s00_dom_libre_seed42, la fila suelo "sin_imputar"
    da ARI=0.000 (1 cluster) con normalizar=True y ARI=1.000 (3 clusters) con
    normalizar=False. Comparar dos metodos que hayan caido en ramas distintas es
    comparar dos pipelines, no dos metodos.

    El resto del pipeline (scale, PCA, vecinos, grid de resolucion) es identico
    para todos, que es lo que mantiene la comparacion limpia.

    Deja constancia de la rama en b.uns["normalizado"] para que quien llama lo
    pueda volcar al CSV sin recalcularlo."""
    b = a.copy()
    auto = es_matriz_de_conteo(a)
    if normalizar is None:
        normalizar = bool(cfg["normalizar"]) and auto
    else:
        normalizar = bool(normalizar)
        if normalizar and not auto:
            # normalize_total/log1p PRESUPONEN conteos no negativos: con
            # negativos log1p(x<-1) da NaN y el PCA revienta. Antes esto se
            # esquivaba cambiando de rama en silencio; ahora se avisa y se hace
            # lo que se ha pedido, que es lo que permite detectarlo.
            X = np.asarray(a.X, dtype=float)
            warnings.warn(
                "construir_grafo: normalizar=True forzado sobre una matriz que "
                "NO es de conteo no negativo (min=%.4f, %.2f%% de negativos). "
                "normalize_total/log1p sobre negativos produce NaN."
                % (np.nanmin(X), 100.0 * float((X < 0).mean())),
                RuntimeWarning, stacklevel=2)
    b.uns["normalizado"] = bool(normalizar)
    if normalizar:
        sc.pp.normalize_total(b, target_sum=1e4)
        sc.pp.log1p(b)
    sc.pp.scale(b, max_value=10)
    n_pcs = min(cfg["n_pcs"], b.n_vars - 1, b.n_obs - 1)
    sc.tl.pca(b, n_comps=n_pcs, random_state=cfg["seed"])
    sc.pp.neighbors(b, n_neighbors=cfg["n_neighbors"], n_pcs=n_pcs,
                    random_state=cfg["seed"])
    return b


def elegir_resolucion(b, n_target, cfg):
    """Corre Leiden en cada resolucion del grid pequeno y devuelve la que da
    n_clusters mas cercano a n_target (desempate -> resolucion mas baja, que da
    la particion mas simple). Devuelve (resolution, labels, n_clusters, barrido)."""
    mejor = None  # (dist_a_target, resolution, labels, n_clusters)
    barrido = []
    for res in cfg["res_grid"]:
        sc.tl.leiden(b, resolution=res, random_state=cfg["seed"], key_added="leiden")
        labels = b.obs["leiden"].astype(int).values
        nclust = len(np.unique(labels))
        barrido.append((res, nclust))
        dist = abs(nclust - n_target)
        if mejor is None or (dist, res) < (mejor[0], mejor[1]):
            mejor = (dist, res, labels, nclust)
    return mejor[1], mejor[2], mejor[3], barrido


def evaluar(dom_real, leiden):
    return dict(
        ARI  = adjusted_rand_score(dom_real, leiden),
        NMI  = normalized_mutual_info_score(dom_real, leiden),
        AMI  = adjusted_mutual_info_score(dom_real, leiden),
        HOMO = homogeneity_score(dom_real, leiden),
        n_clusters = int(len(np.unique(leiden))),
    )


def fig_clustering(coords, g_size, dom_real, leiden, nombre, metricas, outdir,
                   matriz="observado"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].imshow(dom_real.reshape(g_size, g_size), origin="lower", cmap="Set2")
    axes[0].set_title("dominio real (ground truth)")
    axes[1].imshow(leiden.reshape(g_size, g_size), origin="lower", cmap="tab10")
    axes[1].set_title(f"clusters Leiden (n={metricas['n_clusters']})")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(
        f"{nombre} [{matriz}]   ARI={metricas['ARI']:.3f}  NMI={metricas['NMI']:.3f}  "
        f"AMI={metricas['AMI']:.3f}  HOMO={metricas['HOMO']:.3f}",
        fontsize=12)
    fig.tight_layout()
    # figura DE UN DATASET -> subcarpeta de su base
    p = os.path.join(dir_base(nombre, outdir, crear=True),
                     f"fig_clustering_{nombre}_{matriz}.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def main(datasets=None):
    warnings.filterwarnings("ignore")  # silencia FutureWarnings de scanpy/leiden
    datasets = datasets or DATASETS
    if not datasets:
        raise SystemExit(f"No hay datasets con dominios en {OUTDIR}/. "
                         f"Corre antes: python gen_datasets3.py")
    filas = []
    for nombre in datasets:
        for etiqueta_matriz, capa in MATRICES:
            a = cargar(nombre, capa=capa)
            dom_real = a.obs["dominio"].astype(int).values
            g_size = int(round(np.sqrt(a.n_obs)))
            coords = a.obsm["spatial"]

            n_target = int(len(np.unique(dom_real)))
            b = construir_grafo(a, CFG)
            res, leiden, nclust, barrido = elegir_resolucion(b, n_target, CFG)

            m = evaluar(dom_real, leiden)
            # Se registra SI se normalizo: no todas las filas pasan por el
            # mismo preprocesado (las verosimilitudes continuas no se
            # log-normalizan, ver es_matriz_de_conteo) y comparar ARI entre
            # filas con distinto valor aqui es comparar dos pipelines, no dos
            # datasets. Medido: con normalizar=False el colapso de la fila
            # "observado" desaparece en TODOS los datasets (gp_crudo 0.000 ->
            # 1.000, nb_s00 0.000 -> 1.000), asi que la columna no es un
            # detalle: es la variable que explica la mitad de la tabla.
            m.update(dataset=nombre,
                     matriz=etiqueta_matriz,
                     normalizado=bool(CFG["normalizar"] and es_matriz_de_conteo(a)),
                     n_dominios=n_target,
                     sparsity_pct=round(100 * float((np.asarray(a.X) == 0).mean()), 2),
                     resolution=res)
            filas.append(m)

            p = fig_clustering(coords, g_size, dom_real, leiden, nombre, m, OUTDIR,
                               matriz=etiqueta_matriz)
            barrido_str = " ".join(f"r{r}->{n}c" for r, n in barrido)
            print(f"[{nombre:>18} | {etiqueta_matriz:>11}]  norm={str(m['normalizado']):>5}  "
                  f"sparsity={m['sparsity_pct']:>5.2f}%  "
                  f"n_dom={m['n_dominios']}  res_elegida={res}  n_clust={m['n_clusters']:>2}  "
                  f"ARI={m['ARI']:.3f}  NMI={m['NMI']:.3f}  AMI={m['AMI']:.3f}  "
                  f"HOMO={m['HOMO']:.3f}\n{'':>34}barrido: {barrido_str}  -> {os.path.join(base_de(nombre), os.path.basename(p))}")

    cols = ["dataset", "matriz", "normalizado", "sparsity_pct", "n_dominios",
            "n_clusters", "ARI", "NMI", "AMI", "HOMO", "resolution"]
    df = pd.DataFrame(filas)[cols].round(4)
    # tabla GENERICA (compara todas las bases) -> raiz de datasets_gen3/
    csv = os.path.join(OUTDIR, "metricas_clustering.csv")
    df.to_csv(csv, index=False)
    print(f"\nTabla -> {csv}")
    print(df.to_string(index=False))

    # lectura rapida del efecto de la verosimilitud y de la sparsity: media de
    # los dos modos (libre/atada) por dataset base. El regex de antes solo
    # capturaba nb_s\d+ y dejaba fuera gp_crudo/gaussian/poisson, que ahora
    # tambien tienen dominios.
    df["base"] = df["dataset"].str.extract(r"^(.+?)_dom_")
    resumen = df.groupby(["base", "matriz"])[["ARI", "NMI", "AMI", "HOMO"]].mean().round(3)
    print("\nMedia por dataset base (efecto de la verosimilitud y de la sparsity "
          "en la recuperacion de dominios):")
    print(resumen.to_string())

    # margen que deja el held-out: cuanto ARI se pierde por enmascarar y, por
    # tanto, cuanto puede como maximo recuperar un imputador en cada dataset
    piv = df.pivot(index="dataset", columns="matriz", values="ARI")
    if {"observado", "sin_heldout"}.issubset(piv.columns):
        piv["margen_imputacion"] = (piv["sin_heldout"] - piv["observado"]).round(4)
        print("\nMargen de ARI que abre el held-out (techo sin_heldout - suelo observado):")
        print(piv.to_string())


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", nargs="+", default=None,
                    help="datasets a evaluar; por defecto TODOS los "
                         "*_dom_* de datasets_gen3/ (sin las replicas _seedN)")
    args = ap.parse_args()
    main(datasets=args.datasets)
