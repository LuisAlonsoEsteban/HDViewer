r"""
moran.py
========
Calculos de AUTOCORRELACION ESPACIAL (Moran's I), compartidos por los dos
scripts que los usan:

  analyze_datasets.py : Moran GLOBAL de cada dataset, con p-value por
                        permutacion (TABLA 6) -- cuanta estructura espacial
                        SOBREVIVE al ruido de cada verosimilitud.
  gen_datasets3.py    : Moran INTERNO de cada dominio espacial
                        (fig_moran_por_dominio) -- coherencia dentro de cada
                        region, no del campo global.

Definicion (una sola, usada por todo el modulo):
    I = (n/S0) * sum_ij w_ij (x_i - xbar)(x_j - xbar) / sum_i (x_i - xbar)^2

con W = pesos k-NN row-standardized (cada fila suma 1 -> S0 = n, el estandar
en Moran's I). Toda funcion de aqui construye W con build_knn_weights, asi que
un k dado siempre significa lo mismo en los dos scripts.

NOTA sobre el k-NN: se calcula con argsort sobre la matriz densa de distancias,
no con scipy.spatial.cKDTree. Con 400 spots la matriz densa es trivial (400x400)
y evita que el desempate entre vecinos EQUIDISTANTES (frecuente en una rejilla
regular, sobre todo en los spots del borde) dependa del orden en que el arbol
los devuelva -- lo que cambiaria los valores publicados en la TABLA 6.
"""

import numpy as np


# --------------------------------------------------------------------------- #
# PESOS ESPACIALES                                                             #
# --------------------------------------------------------------------------- #
def build_knn_weights(coords, k=8):
    """Matriz de pesos espaciales W (n,n): binaria, k vecinos mas cercanos,
    row-standardized (cada fila suma 1, el estandar en Moran's I)."""
    n = coords.shape[0]
    d2 = ((coords[:, None, :] - coords[None, :, :])**2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    orden = np.argsort(d2, axis=1)[:, :k]
    W = np.zeros((n, n))
    for i in range(n):
        W[i, orden[i]] = 1.0 / k
    return W


# --------------------------------------------------------------------------- #
# MORAN'S I                                                                    #
# --------------------------------------------------------------------------- #
def morans_i_batch(X, W):
    """
    Moran's I para cada columna (gen) de X = (n_spots, n_genes), vectorizado:
        I = (n/S0) * sum_ij w_ij*(x_i-xbar)(x_j-xbar) / sum_i (x_i-xbar)^2
    Devuelve un array (n_genes,) con un I por gen.
    """
    n = X.shape[0]
    Z = X - X.mean(axis=0, keepdims=True)
    WZ = W @ Z
    num = np.sum(Z * WZ, axis=0)
    den = np.sum(Z**2, axis=0)
    S0 = W.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        I = (n / S0) * (num / den)
    return I


def moran_i(valores, coords, k=8):
    """Moran's I de UN solo vector (n_spots,) sobre sus coordenadas.
    Atajo sobre morans_i_batch para cuando no hay varios genes que procesar."""
    valores = np.asarray(valores, dtype=float).reshape(-1, 1)
    W = build_knn_weights(np.asarray(coords, dtype=float), k=k)
    return float(morans_i_batch(valores, W)[0])


def moran_por_dominio(valores, dominio_labels, coords, k=8):
    """Moran's I POR SEPARADO dentro de cada dominio (coherencia espacial
    interna, no del campo global). Devuelve {dominio: moran_medio_genes}.

    valores: (n_spots, n_genes) -- lambda o counts, segun lo que se quiera medir.
    """
    valores = np.asarray(valores, dtype=float)
    coords = np.asarray(coords, dtype=float)
    resultados = {}
    for d in np.unique(dominio_labels):
        m = dominio_labels == d
        # se necesita un minimo de puntos para que el Moran k-NN tenga sentido;
        # con la particion por rechazo esto ya no deberia dispararse, pero se
        # mantiene como red de seguridad. k se reduce a los vecinos disponibles.
        if m.sum() < 5:
            resultados[int(d)] = np.nan
            continue
        # W se construye UNA vez por dominio y se reutiliza para los n_genes
        # genes (morans_i_batch vectoriza sobre columnas)
        W = build_knn_weights(coords[m], k=min(k, int(m.sum()) - 1))
        resultados[int(d)] = float(np.nanmean(morans_i_batch(valores[m], W)))
    return resultados


# --------------------------------------------------------------------------- #
# TEST DE SIGNIFICACION                                                        #
# --------------------------------------------------------------------------- #
def test_moran_permutacion(X, W, n_perm=200, seed=0):
    """
    Test de PERMUTACION para Moran's I (distribution-free: no asume
    normalidad, valido tambien para Poisson/NB).

    H0: no hay autocorrelacion espacial (los valores estan repartidos al azar
    entre spots, sin relacion con su posicion).
    H1: hay autocorrelacion espacial POSITIVA (spots vecinos se parecen) --
    test de una cola, que es lo que predice el kernel RBF.

    Se permuta el ORDEN de los spots (rompe la posicion, mantiene los
    valores) y se recalcula Moran's I cada vez -> distribucion nula empirica.
    """
    rng = np.random.default_rng(seed)
    I_obs = float(np.nanmean(morans_i_batch(X, W)))
    n = X.shape[0]
    nulos = np.empty(n_perm)
    for b in range(n_perm):
        perm = rng.permutation(n)
        nulos[b] = float(np.nanmean(morans_i_batch(X[perm], W)))
    p = float((np.sum(nulos >= I_obs) + 1) / (n_perm + 1))   # una cola, +1 correccion continuidad
    return I_obs, p
