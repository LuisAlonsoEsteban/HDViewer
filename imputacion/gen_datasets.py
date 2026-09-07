r"""
gen_datasets.py
================
SCRIPT 1 de 2: genera los datasets sinteticos + figuras. El analisis
estadistico esta en analyze_datasets.py.

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

Uso:
    python gen_datasets.py
"""

import numpy as np
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import os

from rutas import DATASETS

# --------------------------------------------------------------------------- #
# CONFIG                                                                       #
# --------------------------------------------------------------------------- #
CFG = dict(
    grid_size    = 20,      # 20x20 = 400 spots -> el GP se muestrea en 400 dims
    n_genes      = 100,     # cada gen = 1 realizacion independiente del GP
    gp_mean      = 50.0,    # funcion media mu 
    gp_sigma     = 8.0,     # amplitud del RBF (s_f del lecture = sigma^2)
    length_scale = 3.0,     # length-scale ell del RBF
    jitter       = 1e-6,    # estabilidad numerica de Cholesky
    mask_frac    = 0.10,    # 10% de ground truth 
    lambda_floor = 1e-3,    # lambda > 0 para ser tasa de Poisson/NB
    gauss_tau    = None,    # sd del ruido gaussiano. None -> tau^2 = mu
    incluir_gp_crudo = True,   # guardar tambien X = lambda sin ruido (referencia)
    seed         = 42,
    outdir       = DATASETS,
    target_sparsities = [0.0, 0.10, 0.50],
)


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


# --------------------------------------------------------------------------- #
# VEROSIMILITUDES                                                              #
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
    """
    if modelo == "gp_crudo":
        return lam.copy()
    if modelo == "gaussian":
        return lam + rng.normal(0.0, tau, size=lam.shape)
    if modelo == "poisson":
        return rng.poisson(lam).astype(float)
    if modelo == "nbinom":
        return rng.negative_binomial(r, r / (r + lam)).astype(float)
    raise ValueError(f"modelo desconocido: {modelo!r}")


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
def save(nombre, cfg, coords, counts, obs, mask, lam, rows, cols, gt, extra=None):
    adata = ad.AnnData(X=obs.astype(np.float32))
    adata.obs_names = [f"spot_{i}" for i in range(coords.shape[0])]
    adata.var_names = [f"gene_{g}" for g in range(cfg["n_genes"])]
    adata.obsm["spatial"] = coords
    adata.layers["ground_truth_full"] = counts.astype(np.float32)
    adata.layers["lambda_field"] = lam.astype(np.float32)
    adata.layers["mask"] = mask.astype(np.uint8)
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


# --------------------------------------------------------------------------- #
# FIGURAS                                                                      #
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
    F = sample_gp(K, cfg["gp_mean"], cfg["n_genes"], cfg["jitter"], rng) #ruido gaussiano
    lam = np.clip(F, cfg["lambda_floor"], None) #solo valores positivos
    print(f"[campo] lambda: min={lam.min():.2f} max={lam.max():.2f} "
          f"media={lam.mean():.2f}  clippeados={int((F < cfg['lambda_floor']).sum())}\n")

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
        h5 = save(modelo, cfg, coords, counts, obs, mask, lam, rows, cols, gt, extra)
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
                  {"nb_dispersion_r": float(r), "target_sparsity": float(target)})
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
    print("Siguiente paso:  python analyze_datasets.py")


if __name__ == "__main__":
    main(CFG)
