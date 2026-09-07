r"""
rutas.py
========
UNICA fuente de verdad de las rutas del proyecto.

Todo se deriva de la ubicacion de ESTE fichero (imputacion/rutas.py), no del
directorio desde el que se lance el script. Asi da igual si corres desde
imputacion/, desde synthetic_out/ o desde methods/: las rutas son siempre las
mismas.

Antes habia cuatro convenciones distintas conviviendo para las MISMAS dos
carpetas -- relativas ("datasets", "datasets_gen3", que solo funcionaban con
el CWD en synthetic_out/) y absolutas escritas a mano. Las absolutas de
methods/ apuntaban a imputacion/datasets, que dejo de existir al mover la
carpeta dentro de synthetic_out/, y por eso esos scripts quedaron rotos.
Centralizar aqui evita que vuelva a pasar: si la carpeta se mueve, se cambia
UNA linea.

Uso:
    from rutas import DATASETS as INDIR          # scripts de imputacion/
    from rutas import GEN3 as OUTDIR

    # desde methods/ (subdirectorio: hay que anadir el padre al path antes)
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from rutas import DATASETS as INDIR
"""

import os

# imputacion/
BASE = os.path.dirname(os.path.abspath(__file__))

# imputacion/synthetic_out/
SALIDA = os.path.join(BASE, "synthetic_out")

# datasets de gen_datasets.py / gen_datasets2.py (verosimilitudes: gp_crudo,
# gaussian, poisson, nb_s00, nb_s10, nb_s50) + sus analisis y benchmarks
DATASETS = os.path.join(SALIDA, "datasets")

# datasets de gen_datasets3.py (dominios espaciales: *_dom_libre / *_dom_atada)
# + el clustering de cluster_domains.py
GEN3 = os.path.join(SALIDA, "datasets_gen3")

# imputacion/methods/ -- los impute_*.py. Se expone para que los scripts de
# imputacion/ puedan importarlos (p.ej. evaluar_clustering_knn.py necesita
# impute_knn_spatial) sin depender de PYTHONPATH
METHODS = os.path.join(BASE, "methods")


# --------------------------------------------------------------------------- #
# LAYOUT DE datasets_gen3/: UNA SUBCARPETA POR BASE                            #
# --------------------------------------------------------------------------- #
# datasets_gen3/
#   |- <base>/                     <- todo lo que es DE UNA SOLA BASE
#   |    |- dataset_<base>_dom_<modo>.h5ad
#   |    |- ground_truth_<base>_dom_<modo>.npz
#   |    |- fig_dominios_mapa_<base>.png, fig_moran_por_dominio_<base>.png,
#   |       fig_correlacion_por_variante_<base>.png, fig_clustering_*.png
#   |- metricas_estadisticas.csv   <- lo GENERICO (compara TODAS las bases entre
#   |- metricas_clustering.csv        si) se queda en la raiz
#   |- benchmark_*.csv
#   |- fig_descomposicion_varianza.png, fig_media_varianza.png,
#      fig_dispersion_vs_sparsity.png, fig_tipos_de_cero.png,
#      fig_matriz_confusion.png
#
# Ningun script construye estas rutas a mano: todos pasan por las funciones de
# abajo, para que mover la carpeta o cambiar el layout siga siendo un cambio de
# un solo sitio (que es para lo que existe este fichero).

def base_de(nombre):
    """Base de la que sale un dataset de gen3.

        nb_s10_dom_libre        -> nb_s10
        nb_s10_dom_libre_seed7  -> nb_s10   (las replicas van con su base)
        poisson                 -> poisson  (dataset de gen2, sin capa de dominios)
    """
    return nombre.split("_dom_")[0]


def dir_base(nombre, raiz=None, crear=False):
    """Subcarpeta que le toca a ese dataset: <raiz>/<base>/"""
    raiz = raiz if raiz is not None else GEN3
    d = os.path.join(raiz, base_de(nombre))
    if crear:
        os.makedirs(d, exist_ok=True)
    return d


def _resolver(raiz, fichero, nombre):
    """Busca el fichero en <raiz>/<base>/ y, si no esta, en la raiz plana.

    El fallback NO es transitorio: datasets/ (gen_datasets2) no tiene
    subcarpetas -- alli cada 'base' ES el dataset y una subcarpeta por dataset
    no aporta nada -- y los mismos lectores (analyze_datasets.py) sirven para
    las dos carpetas. De paso, sigue encontrando ficheros de gen3 que se
    hubieran quedado en el layout plano anterior.
    """
    p = os.path.join(raiz, base_de(nombre), fichero)
    return p if os.path.exists(p) else os.path.join(raiz, fichero)


def ruta_dataset(nombre, raiz=None):
    raiz = raiz if raiz is not None else GEN3
    return _resolver(raiz, f"dataset_{nombre}.h5ad", nombre)


def ruta_ground_truth(nombre, raiz=None):
    raiz = raiz if raiz is not None else GEN3
    return _resolver(raiz, f"ground_truth_{nombre}.npz", nombre)


def datasets_en(raiz=None, patron=""):
    """Nombres de los dataset_*.h5ad que hay en raiz, mirando TAMBIEN dentro de
    las subcarpetas de base. Devuelve nombres (sin 'dataset_' ni '.h5ad'), sin
    duplicados: es lo que usan los scripts para descubrir que hay generado."""
    raiz = raiz if raiz is not None else GEN3
    nombres = set()
    if not os.path.isdir(raiz):
        return []
    for carpeta, _, ficheros in os.walk(raiz):
        for f in ficheros:
            if f.startswith("dataset_") and f.endswith(".h5ad") and patron in f:
                nombres.add(f[len("dataset_"):-len(".h5ad")])
    return sorted(nombres)
