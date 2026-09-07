r"""
gen_datasets.py
================
SCRIPT UNICO: genera los datasets sinteticos + figuras de generacion (campo GP,
kernel, verosimilitudes, mascara) Y ADEMAS las figuras de distribucion
(histogramas vs teorica, ECDF, Q-Q, violin). Antes eran dos scripts
(gen_datasets.py + plot_distributions.py); ahora todo corre en un solo main().

=============================================================================
DISENO: comparar MODELOS DE OBSERVACION (verosimilitudes) sobre el MISMO campo
=============================================================================
Un unico campo latente lambda, generado por el GP, COMPARTIDO por todos los
datasets (mismo seed). Lo unico que cambia es la VEROSIMILITUD que se aplica
encima:

    lambda ~ GP(mu=50, K_RBF)                 <- campo latente, COMUN a todos
      |
      +-- gaussian : X = lambda + eps,  eps ~ N(0, tau^2)   <- verosim. GAUSSIANA
      +-- poisson  : X ~ Poisson(lambda)                     <- verosim. POISSON
      +-- nb_sXX   : X ~ NegBinom(media=lambda, disp=r)      <- verosim. NB
                     con r ajustado por biseccion -> sparsity 0%, 10%, 50%

Sigue el modelo del lecture (pag. 95):
    y = f(x) + eps ,  eps ~ N(0, sigma_n^2)
donde sigma_n es la "noise standard deviation" que aparece como K_XX + sigma^2*I
(lecture, pags. 75-76).

Cada dataset tiene EXACTAMENTE UNA verosimilitud. NUNCA se mete eps Y Poisson
a la vez en el mismo dataset: eso seria solapar los dos ruidos. Asi la
comparacion es limpia: mismo lambda, una sola verosimilitud cada vez.

Matiz que NO es solapamiento sino la estructura del modelo: lambda viene de un
GP, luego es aleatorio. Eso es INEVITABLE y es justo lo que da la estructura
espacial. Es un modelo JERARQUICO de dos niveles:
    NIVEL 1: lambda ~ GP        -> la SENAL espacial (predecible desde vecinos)
    NIVEL 2: X | lambda ~ ...   -> el RUIDO de observacion (impredecible)
La varianza total mezcla ambos:  Var(X) = Var(lambda) + E[Var(X|lambda)].
analyze_datasets.py los separa numericamente.

TAU: por defecto tau^2 = mu, para que el ruido gaussiano tenga la MISMA
varianza media que tendria el de Poisson. Asi la comparacion aisla la FORMA del
ruido (homocedastico vs acoplado a la media), no su cantidad.

--------------------------------------------------------------------------- #
NOTA sobre numpy vs scipy.stats en el muestreo
--------------------------------------------------------------------------- #
El campo GP (sample_gp) se queda con Cholesky manual (L @ Z): es exactamente
la formula del lecture (f = mu + L@Z) y para una matriz 400x400 no gana nada
usar scipy.stats.multivariate_normal (que internamente hace SVD/eigh, no
Cholesky).

Las verosimilitudes de observacion (gaussian, poisson, nbinom) SI pasan a
scipy.stats (norm.rvs, poisson.rvs, nbinom.rvs) con random_state=rng. Motivo:
las figuras de distribucion (mas abajo) ya validan el muestreo contra
scipy.stats.norm.pdf / poisson.pmf / nbinom.pmf. Usar los mismos objetos de
scipy.stats para generar Y para validar evita tener la parametrizacion de la
NB escrita dos veces (una para muestrear, otra para el pmf teorico) con riesgo
de que diverjan. La parametrizacion (n=r, p=r/(r+lam)) es identica a la de
numpy.random.negative_binomial, asi que el cambio no altera el comportamiento
estadistico, solo unifica la fuente de verdad.

Uso:
    python gen_datasets.py
"""

import os
import numpy as np
from scipy.spatial.distance import cdist
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

from rutas import DATASETS

# --------------------------------------------------------------------------- #
# CONFIG                                                                       #
# --------------------------------------------------------------------------- #
CFG = dict(
    grid_size    = 20,
    n_genes      = 100,
    gp_mean      = 50.0,
    gp_sigma     = 8.0,
    length_scale = 3.0,
    jitter       = 1e-6,
    mask_frac    = 0.10,
    lambda_floor = 1e-3,
    gauss_tau    = None,
    incluir_gp_crudo = True,
    seed         = 42,
    outdir       = DATASETS,
    target_sparsities = [0.0, 0.10, 0.50],
    correlacion_genica = True,      # <- sin comillas, con =
    n_factores         = 5,
    genes_por_modulo   = None,
    peso_intra         = 1.0,
    peso_fondo         = 0.15,
)

ORDEN = ["gp_crudo", "gaussian", "poisson", "nb_s00", "nb_s10", "nb_s50"]
COLORES = dict(gp_crudo="#555555", gaussian="seagreen", poisson="steelblue",
               nb_s00="#999999", nb_s10="darkorange", nb_s50="firebrick")


# --------------------------------------------------------------------------- #
# GRID + KERNEL + GP                                                           #
# --------------------------------------------------------------------------- #
def build_grid(n):
    """Entrada 2D del GP: coordenadas (x, y) de los spots."""
    xs, ys = np.meshgrid(np.arange(n), np.arange(n))
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def rbf_kernel(coords, length_scale, sigma):
    """
    Covarianza del GP (lecture, pag. 82):
        k(x, x') = s_f * exp( -1/2 * (x-x')^2 / ell^2 )
    con s_f = sigma^2. Equivale a sklearn.gaussian_process.kernels.RBF
    escalado.
    """
    sq = cdist(coords, coords, metric="sqeuclidean")
    return sigma**2 * np.exp(-sq / (2.0 * length_scale**2))


def sample_gp(K, mean, n_genes, jitter, rng):
    """
    Muestreo del GP: f ~ N(mu, K) via Cholesky K = L L'.
        F = mu + L @ Z,  Z ~ N(0, I)
    Cada COLUMNA de F es una realizacion de 400 dims (= un gen).
    Todas comparten la MISMA covarianza K.
    """
    n = K.shape[0]
    L = np.linalg.cholesky(K + jitter * np.eye(n))
    return mean + L @ rng.standard_normal((n, n_genes))

def sample_spatial_factors(K, n_factores, jitter, rng):
    """
    Genera n_factores campos GP espaciales INDEPENDIENTES (los 'programas
    transcripcionales' o modulos). Cada factor es una realizacion del GP
    con la misma covarianza espacial K. Devuelve (n_spots, n_factores).
    """
    n_spots = K.shape[0]
    L = np.linalg.cholesky(K + jitter * np.eye(n_spots))
    Z = rng.standard_normal((n_spots, n_factores))
    factores = L @ Z          # (n_spots, n_factores), media 0
    return factores


def build_loadings(n_factores, n_genes, genes_por_modulo, peso_intra,
                   peso_fondo, rng):
    """
    Matriz de cargas W (n_factores, n_genes) que asigna genes a modulos.

    - Cada gen se asocia PRINCIPALMENTE a un factor (su modulo), con carga
      alta (peso_intra) -> genes del mismo modulo quedan correlacionados.
    - Ademas recibe cargas pequenas y aleatorias sobre el resto de factores
      (peso_fondo) -> co-expresion parcial realista, no bloques perfectos.

    genes_por_modulo : lista con cuantos genes tiene cada modulo. Debe sumar
                       n_genes. Si None, se reparten equitativamente.
    """
    if genes_por_modulo is None:
        base = n_genes // n_factores
        genes_por_modulo = [base] * n_factores
        genes_por_modulo[-1] += n_genes - sum(genes_por_modulo)  # resto al ultimo
    assert sum(genes_por_modulo) == n_genes, "genes_por_modulo debe sumar n_genes"

    W = peso_fondo * rng.standard_normal((n_factores, n_genes))  # fondo debil
    col = 0
    asignacion = np.zeros(n_genes, dtype=int)  # a que modulo pertenece cada gen
    for f, ngm in enumerate(genes_por_modulo):
        # carga fuerte del factor f sobre sus genes (con algo de variabilidad)
        W[f, col:col + ngm] += peso_intra * (0.5 + rng.random(ngm))
        asignacion[col:col + ngm] = f
        col += ngm
    return W, asignacion


def sample_gp_correlacionado(K, mu, n_genes, n_factores, genes_por_modulo,
                             peso_intra, peso_fondo, jitter, rng):
    """
    Campo latente lambda CON correlacion espacial Y genica.

    lambda[:, g] = mu + sum_f W[f, g] * factor_espacial_f

    Genes que cargan sobre el mismo factor comparten patron espacial ->
    correlacionados entre si. Factores distintos -> patrones espaciales
    independientes -> genes de modulos distintos no correlacionados.

    Devuelve (F, W, asignacion):
      F           : (n_spots, n_genes) campo latente
      W           : (n_factores, n_genes) cargas usadas
      asignacion  : (n_genes,) modulo de cada gen (para validacion posterior)
    """
    factores = sample_spatial_factors(K, n_factores, jitter, rng)  # (n_spots, n_fact)
    W, asignacion = build_loadings(n_factores, n_genes, genes_por_modulo,
                                   peso_intra, peso_fondo, rng)
    F = mu + factores @ W      # (n_spots, n_genes)
    return F, W, asignacion, factores
# --------------------------------------------------------------------------- #
# VEROSIMILITUDES (con scipy.stats, misma parametrizacion que la validacion)   #
# --------------------------------------------------------------------------- #
def nb_zero_prob(lam, r):
    """P(X=0) media de una NB(media=lam, disp=r):  (r/(r+lam))^r"""
    return float(np.mean((r / (r + lam))**r))


def find_r_for_sparsity(lam, target, lo=1e-4, hi=1e7, tol=1e-6, max_iter=200):
    """
    Biseccion en escala log sobre r para alcanzar una sparsity objetivo.
    P(X=0) es monotona DECRECIENTE en r.
    target=0 es inalcanzable exacto (haria falta r -> inf): se devuelve r muy
    grande = limite de Poisson (con mu=50 la sparsity es e^-50 ~ 0).
    """
    if target <= 0.0:
        return hi, nb_zero_prob(lam, hi)
    for _ in range(max_iter):
        mid = np.sqrt(lo * hi)
        p = nb_zero_prob(lam, mid)
        if abs(p - target) < tol:
            return mid, p
        if p > target:
            lo = mid
        else:
            hi = mid
    return mid, nb_zero_prob(lam, mid)


def observe(lam, modelo, rng, tau=None, r=None):
    """
    Aplica UNA verosimilitud al campo lambda. Nunca dos a la vez.

      "gp_crudo": X = lambda                       (sin ruido, referencia)
      "gaussian": X = lambda + N(0, tau^2)         Var(X|lam) = tau^2  CONSTANTE
      "poisson" : X ~ Poisson(lambda)              Var(X|lam) = lambda (crece con la media)
      "nbinom"  : X ~ NegBinom(media=lam, disp=r)  Var(X|lam) = lam + lam^2/r

    Muestreo via scipy.stats (norm/poisson/nbinom .rvs con random_state=rng)
    para compartir parametrizacion con las figuras de distribucion, que
    validan contra .pdf/.pmf de las mismas familias.
    """
    if modelo == "gp_crudo":
        return lam.copy()
    if modelo == "gaussian":
        return stats.norm.rvs(loc=lam, scale=tau, random_state=rng)
    if modelo == "poisson":
        return stats.poisson.rvs(mu=lam, random_state=rng).astype(float)
    if modelo == "nbinom":
        return stats.nbinom.rvs(r, r / (r + lam), random_state=rng).astype(float)
    raise ValueError(f"modelo desconocido: {modelo!r}")

def validar_correlacion_genica(lam, asignacion, outdir):
    """
    Valida y visualiza la estructura de co-expresion genica del campo latente.
    Genera un heatmap de la matriz de correlacion gen-gen (ordenado por modulo,
    deberian verse bloques diagonales) y calcula el ratio de correlacion media
    intra-modulo vs inter-modulo.

    Devuelve (ruta_figura, ratio_intra_inter).
    """
    import matplotlib.pyplot as plt

    orden = np.argsort(asignacion)         # agrupa genes contiguos por modulo
    asign_ord = asignacion[orden]
    C = np.corrcoef(lam[:, orden].T)       # (n_genes, n_genes)

    # correlacion media intra vs inter modulo (sobre la diagonal superior)
    n = len(asignacion)
    iu = np.triu_indices(n, k=1)
    corr_pairs = np.corrcoef(lam.T)[iu]
    mismo_modulo = (asignacion[:, None] == asignacion[None, :])[iu]
    intra = corr_pairs[mismo_modulo]
    inter = corr_pairs[~mismo_modulo]
    intra_m = float(np.mean(intra)) if len(intra) else float("nan")
    inter_m = float(np.mean(inter)) if len(inter) else float("nan")
    ratio = intra_m / abs(inter_m) if inter_m != 0 else float("inf")

    print(f"[corr]  intra-modulo media={intra_m:.3f}  "
          f"inter-modulo media={inter_m:.3f}  ratio={ratio:.2f}  "
          f"(quieres intra >> inter)")

    # heatmap con lineas separando modulos
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    limites = np.where(np.diff(asign_ord) != 0)[0] + 0.5
    for lb in limites:
        ax.axhline(lb, color="black", lw=0.6)
        ax.axvline(lb, color="black", lw=0.6)
    ax.set_title(f"Correlacion gen-gen (ordenado por modulo)\n"
                 f"intra={intra_m:.2f} inter={inter_m:.2f} ratio={ratio:.1f}")
    ax.set_xlabel("genes"); ax.set_ylabel("genes")
    plt.colorbar(im, label="correlacion de Pearson")
    p = os.path.join(outdir, "correlacion_genica.png")
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close()
    return p, ratio

# --------------------------------------------------------------------------- #
# ENMASCARADO -> GROUND TRUTH                                                  #
# --------------------------------------------------------------------------- #
def mask_ground_truth(counts, mask_frac, rng):
    """
    Enmascara mask_frac de los genes de cada spot y guarda los valores reales.
    OJO: una entrada enmascarada puede valer 0 de por si (cero biologico, solo
    en los datasets de NB). Se registra igual: el metodo debe aprender a
    devolver 0 ahi tambien.
    """
    n_spots, n_genes = counts.shape
    k = int(round(mask_frac * n_genes))
    mask = np.zeros((n_spots, n_genes), dtype=bool)
    for i in range(n_spots):
        mask[i, rng.choice(n_genes, size=k, replace=False)] = True
    obs = counts.copy()
    obs[mask] = 0.0
    rows, cols = np.where(mask)
    return obs, mask, rows, cols, counts[mask].copy()


# --------------------------------------------------------------------------- #
# GUARDADO                                                                     #
# --------------------------------------------------------------------------- #
def save(nombre, cfg, coords, counts, obs, mask, lam, rows, cols, gt,
         asignacion=None, W=None, factores=None, extra=None):
    adata = ad.AnnData(X=obs.astype(np.float32))
    adata.obs_names = [f"spot_{i}" for i in range(coords.shape[0])]
    adata.var_names = [f"gene_{g}" for g in range(cfg["n_genes"])]
    adata.obsm["spatial"] = coords
    adata.layers["ground_truth_full"] = counts.astype(np.float32)
    adata.layers["lambda_field"] = lam.astype(np.float32)
    adata.layers["mask"] = mask.astype(np.uint8)

    # --- trazabilidad de la estructura de co-expresion (si aplica) ---
    if asignacion is not None:
        adata.var["modulo"] = asignacion
    if W is not None:
        # W: (n_factores, n_genes). La guardamos en varm (indexado por gen)
        # transponiendola -> (n_genes, n_factores), un vector de cargas por gen.
        adata.varm["loadings"] = W.T.astype(np.float32)
    if factores is not None:
        # factores: (n_spots, n_factores). Los campos GP espaciales latentes
        # de cada modulo, indexados por spot -> obsm.
        adata.obsm["factores_espaciales"] = factores.astype(np.float32)

    params = {k: v for k, v in cfg.items()
              if k not in ("outdir", "target_sparsities", "gauss_tau",
                           "incluir_gp_crudo")}
    params["modelo"] = nombre
    if extra:
        params.update(extra)
    adata.uns["generation_params"] = params
    h5 = os.path.join(cfg["outdir"], f"dataset_{nombre}.h5ad")
    adata.write_h5ad(h5)
    np.savez_compressed(os.path.join(cfg["outdir"], f"ground_truth_{nombre}.npz"),
                        spot_idx=rows, gene_idx=cols, true_value=gt)
    return h5


def cargar(nombre, indir):
    """Recarga un dataset guardado (usado por las figuras de distribucion)."""
    a = ad.read_h5ad(os.path.join(indir, f"dataset_{nombre}.h5ad"))
    p = dict(a.uns["generation_params"])
    counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
    lam = np.asarray(a.layers["lambda_field"], dtype=float)
    return counts, lam, p


# --------------------------------------------------------------------------- #
# FIGURAS DE GENERACION                                                        #
# --------------------------------------------------------------------------- #
def fig_campo_gp(lam, g_size, outdir):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for g, ax in enumerate(axes):
        im = ax.imshow(lam[:, g].reshape(g_size, g_size), origin="lower", cmap="viridis")
        ax.set_title(f"gen {g}"); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("Campo latente lambda del GP -- COMUN a todos los datasets\n"
                 "cada gen = 1 realizacion independiente del mismo proceso", fontsize=12)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_campo_gp.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_kernel(coords, F, ls, sigma, outdir):
    d = cdist(coords, coords).ravel()
    Fc = F - F.mean(axis=0, keepdims=True)
    emp = (Fc @ Fc.T / F.shape[1]).ravel()
    bins = np.linspace(0, d.max(), 30)
    idx = np.digitize(d, bins)
    xs, ys = [], []
    for b in range(1, len(bins)):
        m = idx == b
        if m.sum():
            xs.append(d[m].mean()); ys.append(emp[m].mean())
    rr = np.linspace(0, d.max(), 200)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(xs, ys, s=22, label="covarianza empirica")
    ax.plot(rr, sigma**2 * np.exp(-rr**2 / (2 * ls**2)), "r-", lw=2, label="RBF teorico")
    ax.set_xlabel("distancia entre spots"); ax.set_ylabel("covarianza")
    ax.set_title("Validacion del kernel: la covarianza decae como el RBF")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_kernel.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_verosimilitudes(datasets, g_size, outdir):
    """Mismo campo lambda, distintas verosimilitudes."""
    n = len(datasets)
    fig, axes = plt.subplots(2, n, figsize=(3.5 * n, 7.5))
    g = 0
    for j, (nombre, counts) in enumerate(datasets.items()):
        ax = axes[0, j]
        im = ax.imshow(counts[:, g].reshape(g_size, g_size), origin="lower", cmap="viridis")
        ax.set_title(nombre, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046)
        ax = axes[1, j]
        v = counts.ravel()
        ax.hist(v[v < np.quantile(v, 0.995)], bins=50, color="steelblue")
        ax.set_yscale("log")
        neg = 100 * float((counts < 0).mean())
        sp = 100 * float((counts == 0).mean())
        ax.set_title(f"sparsity={sp:.1f}%   negativos={neg:.2f}%", fontsize=9)
        ax.set_xlabel("valor")
        if neg > 0:
            ax.axvline(0, color="red", ls="--", lw=1.5)
    axes[0, 0].set_ylabel("mapa espacial (gen 0)")
    axes[1, 0].set_ylabel("histograma (log)")
    fig.suptitle("MISMO campo lambda, distintas VEROSIMILITUDES\n"
                 "(linea roja = cero; a su izquierda hay valores negativos, imposibles "
                 "en datos de conteo)", fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_verosimilitudes.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_mascara(obs, mask, counts, g_size, outdir, nombre):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    g = 0
    im = axes[0].imshow(counts[:, g].reshape(g_size, g_size), origin="lower", cmap="viridis")
    axes[0].set_title(f"gen {g}: verdad completa")
    fig.colorbar(im, ax=axes[0], fraction=0.046)
    im = axes[1].imshow(obs[:, g].reshape(g_size, g_size), origin="lower", cmap="viridis")
    axes[1].set_title(f"gen {g}: observado (con huecos)")
    fig.colorbar(im, ax=axes[1], fraction=0.046)
    tipo = np.zeros(counts.shape[0])
    tipo[(counts[:, g] == 0) & (~mask[:, g])] = 1
    tipo[mask[:, g]] = 2
    im = axes[2].imshow(tipo.reshape(g_size, g_size), origin="lower",
                        cmap="Set1", vmin=0, vmax=2)
    axes[2].set_title("0=dato   1=cero BIOLOGICO   2=cero de MASCARA", fontsize=9)
    fig.colorbar(im, ax=axes[2], fraction=0.046)
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"[{nombre}] Los dos tipos de cero: recuperar los de mascara, "
                 f"NO tocar los biologicos", fontsize=12)
    fig.tight_layout()
    p = os.path.join(outdir, f"fig_mascara_{nombre}.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


# --------------------------------------------------------------------------- #
# FIGURAS DE DISTRIBUCION (antes plot_distributions.py)                       #
# --------------------------------------------------------------------------- #
def fig_dist_histogramas(indir):
    disponibles = [n for n in ORDEN if os.path.exists(os.path.join(indir, f"dataset_{n}.h5ad"))]
    n = len(disponibles)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.6 * nrow))
    axes = np.array(axes).ravel()

    for ax, nombre in zip(axes, disponibles):
        counts, lam, p = cargar(nombre, indir)
        vals = counts.ravel()
        mu = float(lam.mean())
        modelo = str(p.get("modelo", nombre))

        vmax = np.quantile(vals, 0.995)
        ax.hist(vals[vals <= vmax], bins=50, density=True, alpha=0.6,
                color=COLORES.get(nombre, "gray"), label="empirico")

        if modelo in ("gp_crudo", "gaussian"):
            xs = np.linspace(vals.min(), vmax, 200)
            tau = p.get("gauss_tau", np.sqrt(mu))
            if modelo == "gaussian":
                ax.plot(xs, stats.norm.pdf(xs, mu, tau), "r-", lw=2,
                        label=f"Normal(mu={mu:.0f}, tau={tau:.1f})")
            else:
                ax.plot([], [], " ", label="(campo lambda, sin ruido)")
        else:
            ks = np.arange(0, int(vmax) + 1)
            if modelo == "poisson" or (modelo.startswith("nb")
                                       and p.get("nb_dispersion_r", 0) > 1e5):
                pmf = stats.poisson.pmf(ks, mu)
                lab = f"Poisson(lambda~{mu:.0f})"
            else:
                r = p["nb_dispersion_r"]
                pmf = stats.nbinom.pmf(ks, r, r / (r + mu))
                lab = f"NB(mu~{mu:.0f}, r={r:.3g})"
            ax.plot(ks, pmf, "r-", lw=1.5, label=lab)

        sp = 100 * float((vals == 0).mean())
        ax.set_title(f"{nombre}   (sparsity={sp:.1f}%)", fontsize=10)
        ax.set_xlabel("valor"); ax.set_ylabel("densidad")
        ax.legend(fontsize=8)

    for ax in axes[len(disponibles):]:
        ax.axis("off")
    fig.suptitle("Distribucion de cada dataset: histograma empirico vs distribucion "
                 "teorica (scipy)\nsi coinciden, el muestreo esta bien hecho", fontsize=13)
    fig.tight_layout()
    path = os.path.join(indir, "fig_dist_histogramas.png")
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return path


def fig_dist_ecdf(indir):
    disponibles = [n for n in ORDEN if os.path.exists(os.path.join(indir, f"dataset_{n}.h5ad"))]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for nombre in disponibles:
        counts, _, _ = cargar(nombre, indir)
        vals = np.sort(counts.ravel())
        ecdf = np.arange(1, len(vals) + 1) / len(vals)
        ax1.plot(vals, ecdf, lw=1.8, color=COLORES.get(nombre), label=nombre)
        m = vals <= 100
        ax2.plot(vals[m], ecdf[m], lw=1.8, color=COLORES.get(nombre), label=nombre)

    for ax, tit in [(ax1, "ECDF completa (escala log en X)"),
                    (ax2, "ECDF zoom (0-100): se ven los ceros y la masa central")]:
        ax.set_xlabel("valor"); ax.set_ylabel("F(x) = P(X <= x)")
        ax.set_title(tit, fontsize=10); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax1.set_xscale("symlog")

    fig.suptitle("Funcion de distribucion acumulada empirica (ECDF)\n"
                 "el salto en x=0 = sparsity; la cola derecha = sobredispersion", fontsize=13)
    fig.tight_layout()
    path = os.path.join(indir, "fig_dist_ecdf.png")
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return path


def fig_dist_qq(indir):
    objetivo = [n for n in ["poisson", "nb_s10", "nb_s50"]
                if os.path.exists(os.path.join(indir, f"dataset_{n}.h5ad"))]
    fig, axes = plt.subplots(1, len(objetivo), figsize=(5 * len(objetivo), 4.6))
    if len(objetivo) == 1:
        axes = [axes]

    for ax, nombre in zip(axes, objetivo):
        counts, lam, p = cargar(nombre, indir)
        vals = counts.ravel()
        mu = float(lam.mean())
        modelo = str(p.get("modelo", nombre))

        probs = (np.arange(1, len(vals) + 1) - 0.5) / len(vals)
        if modelo.startswith("nb"):
            r = p["nb_dispersion_r"]
            q_teo = stats.nbinom.ppf(probs, r, r / (r + mu))
        else:
            q_teo = stats.poisson.ppf(probs, mu)
        q_emp = np.sort(vals)

        ax.scatter(q_teo, q_emp, s=4, alpha=0.3, color=COLORES.get(nombre))
        lim = [0, max(np.nanmax(q_teo[np.isfinite(q_teo)]), q_emp.max())]
        ax.plot(lim, lim, "r--", lw=1.5, label="diagonal (ajuste perfecto)")
        ax.set_xlabel("cuantiles teoricos"); ax.set_ylabel("cuantiles empiricos")
        ax.set_title(f"{nombre}", fontsize=11)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.suptitle("Q-Q plot: cuantiles empiricos vs teoricos\n"
                 "si el muestreo es correcto, los puntos caen sobre la diagonal", fontsize=13)
    fig.tight_layout()
    path = os.path.join(indir, "fig_dist_qq_poisson.png")
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return path


def fig_dist_violin(indir):
    disponibles = [n for n in ORDEN if os.path.exists(os.path.join(indir, f"dataset_{n}.h5ad"))]
    datos = []
    for nombre in disponibles:
        counts, _, _ = cargar(nombre, indir)
        v = counts.ravel()
        datos.append(v[v <= np.quantile(v, 0.99)])

    fig, ax = plt.subplots(figsize=(11, 5.5))
    parts = ax.violinplot(datos, showmedians=True, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(COLORES.get(disponibles[i], "gray"))
        pc.set_alpha(0.65)
    ax.set_xticks(range(1, len(disponibles) + 1))
    ax.set_xticklabels(disponibles, rotation=15)
    ax.set_ylabel("valor (recortado al percentil 99 para ver la forma)")
    ax.set_title("Distribucion comparada por dataset (violin, mediana marcada)\n"
                 "mismo campo lambda, distinta verosimilitud -> distinta forma y dispersion",
                 fontsize=12)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    path = os.path.join(indir, "fig_dist_violin.png")
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# MAIN                                                                         #
# --------------------------------------------------------------------------- #
def main(cfg):
    os.makedirs(cfg["outdir"], exist_ok=True)
    coords = build_grid(cfg["grid_size"])
    tau = cfg["gauss_tau"] if cfg["gauss_tau"] else np.sqrt(cfg["gp_mean"])

    print(f"[grid]  {coords.shape[0]} spots ({cfg['grid_size']}x{cfg['grid_size']})"
          f" -> GP muestreado en {coords.shape[0]} dimensiones")
    print(f"[gp]    mu={cfg['gp_mean']}  RBF: sigma={cfg['gp_sigma']} "
          f"(s_f={cfg['gp_sigma']**2:.0f})  ell={cfg['length_scale']}")
    print(f"[genes] {cfg['n_genes']} genes = {cfg['n_genes']} realizaciones del MISMO GP")
    print(f"[mask]  {100*cfg['mask_frac']:.0f}% de ground truth")
    print(f"[tau]   {tau:.3f}  ->  tau^2={tau**2:.1f} = varianza media que tendria "
          f"la Poisson (comparacion justa)\n")

    # --- campo latente lambda: IDENTICO para todos ---
    rng = np.random.default_rng(cfg["seed"])
    K = rbf_kernel(coords, cfg["length_scale"], cfg["gp_sigma"])

    if cfg.get("correlacion_genica", False):
        F, W, asignacion, factores = sample_gp_correlacionado(
            K, cfg["gp_mean"], cfg["n_genes"],
            n_factores=cfg["n_factores"],
            genes_por_modulo=cfg.get("genes_por_modulo", None),
            peso_intra=cfg["peso_intra"],
            peso_fondo=cfg["peso_fondo"],
            jitter=cfg["jitter"], rng=rng,
        )
        print(f"[modulos] {cfg['n_factores']} modulos de co-expresion, "
              f"peso_intra={cfg['peso_intra']}, peso_fondo={cfg['peso_fondo']}")
    else:
        F = sample_gp(K, cfg["gp_mean"], cfg["n_genes"], cfg["jitter"], rng)
        W, asignacion, factores = None, None, None

    lam = np.clip(F, cfg["lambda_floor"], None)
    print(f"[campo] lambda: min={lam.min():.2f} max={lam.max():.2f} "
          f"media={lam.mean():.2f}  clippeados={int((F < cfg['lambda_floor']).sum())}\n")

    # --- validacion de la estructura de co-expresion genica ---
    if cfg.get("correlacion_genica", False) and asignacion is not None:
        p_corr, ratio = validar_correlacion_genica(lam, asignacion, cfg["outdir"])
        print(f"[fig]   {p_corr}\n")
        if ratio < 3.0:
            print(f"[AVISO] ratio intra/inter = {ratio:.2f} < 3.0: los modulos "
                  f"estan poco separados. Considera subir peso_intra o bajar "
                  f"peso_fondo antes de correr los metodos.\n")

    generados = {}
    mascaras = {}

    modelos = ["gaussian", "poisson"]
    if cfg["incluir_gp_crudo"]:
        modelos = ["gp_crudo"] + modelos

    for modelo in modelos:
        r_ = np.random.default_rng(cfg["seed"] + 1)
        counts = observe(lam, modelo, r_, tau=tau)
        obs, mask, rows, cols, gt = mask_ground_truth(counts, cfg["mask_frac"], r_)
        extra = {"gauss_tau": float(tau)} if modelo == "gaussian" else None
        h5 = save(modelo, cfg, coords, counts, obs, mask, lam, rows, cols, gt,
          asignacion=asignacion, W=W, factores=factores, extra=extra)
        generados[modelo] = counts
        mascaras[modelo] = (obs, mask, counts)
        neg = 100 * float((counts < 0).mean())
        print(f"[{modelo:>9}]  sparsity={100*float((counts==0).mean()):>6.3f}%  "
              f"negativos={neg:>6.3f}%  -> {os.path.basename(h5)}")

    print()
    for target in cfg["target_sparsities"]:
        r, p_teo = find_r_for_sparsity(lam, target)
        nombre = f"nb_s{int(round(target*100)):02d}"
        r_ = np.random.default_rng(cfg["seed"] + 1)
        counts = observe(lam, "nbinom", r_, r=r)
        obs, mask, rows, cols, gt = mask_ground_truth(counts, cfg["mask_frac"], r_)
        h5 = save(nombre, cfg, coords, counts, obs, mask, lam, rows, cols, gt,
          asignacion=asignacion, W=W, factores=factores,
          extra={"nb_dispersion_r": float(r), "target_sparsity": float(target)})
        generados[nombre] = counts
        mascaras[nombre] = (obs, mask, counts)
        print(f"[{nombre:>9}]  objetivo={target*100:>4.0f}%  r={r:>11.5g}  "
              f"P(0)teo={p_teo*100:>6.3f}%  real={100*float((counts==0).mean()):>6.3f}%"
              f"  -> {os.path.basename(h5)}")

    print()
    figs = [fig_campo_gp(lam, cfg["grid_size"], cfg["outdir"]),
            fig_kernel(coords, F, cfg["length_scale"], cfg["gp_sigma"], cfg["outdir"]),
            fig_verosimilitudes(generados, cfg["grid_size"], cfg["outdir"])]
    for nm in ["poisson", "nb_s50"]:
        if nm in mascaras:
            figs.append(fig_mascara(*mascaras[nm], cfg["grid_size"], cfg["outdir"], nm))
    for p in figs:
        print(f"[fig]   {p}")

    print(f"\n{len(generados)} datasets en {cfg['outdir']}/")

    # --- figuras de distribucion (antes plot_distributions.py) ---
    print("\n[dist]  generando figuras de distribucion...")
    dist_figs = [fig_dist_histogramas(cfg["outdir"]),
                 fig_dist_ecdf(cfg["outdir"]),
                 fig_dist_qq(cfg["outdir"]),
                 fig_dist_violin(cfg["outdir"])]
    for p in dist_figs:
        print(f"[fig]   {p}")

    print("\nSiguiente paso:  python analyze_datasets.py")


if __name__ == "__main__":
    main(CFG)
