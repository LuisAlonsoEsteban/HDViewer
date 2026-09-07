r"""
figuras.py
==========
TODAS las figuras de validacion, en un solo sitio. Las usan:

  analyze_datasets.py : fig_descomposicion, fig_media_varianza,
                        fig_dispersion_sparsity, fig_ceros, fig_matriz_confusion
  gen_datasets3.py    : fig_dominios_mapa, fig_moran_por_dominio,
                        fig_correlacion_por_variante

Convencion comun a todas: reciben DATOS YA CALCULADOS + un outdir, guardan el
PNG y devuelven su ruta. Ninguna lee ficheros ni calcula estadisticos -- eso
vive en analyze_datasets.py (descomposicion, dispersion, tests) y en moran.py
(autocorrelacion espacial). Asi este modulo no importa nada de los scripts que
lo usan y no hay imports circulares.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# VARIANZA Y VEROSIMILITUD                                                     #
# --------------------------------------------------------------------------- #
def fig_descomposicion(res, outdir):
    nombres = [r["dataset"] for r in res]
    sen = np.array([r["var_senal"] for r in res])
    rui = np.array([r["var_ruido"] for r in res])
    x = np.arange(len(nombres))

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    ax = axes[0]
    ax.bar(x, sen, label="Var SENAL (campo GP)", color="steelblue")
    ax.bar(x, rui, bottom=sen, label="Var RUIDO (verosimilitud)", color="darkorange")
    ax.set_yscale("symlog")
    ax.set_xticks(x); ax.set_xticklabels(nombres, rotation=30, ha="right")
    ax.set_ylabel("varianza (log)")
    ax.set_title("Descomposicion de la varianza (absoluta)")
    ax.legend(fontsize=8)

    ax = axes[1]
    tot = sen + rui
    ps = 100 * sen / np.where(tot > 0, tot, 1)
    pr = 100 * rui / np.where(tot > 0, tot, 1)
    ax.bar(x, ps, label="% SENAL", color="steelblue")
    ax.bar(x, pr, bottom=ps, label="% RUIDO", color="darkorange")
    ax.set_xticks(x); ax.set_xticklabels(nombres, rotation=30, ha="right")
    ax.set_ylabel("% de la varianza total"); ax.set_ylim(0, 100)
    ax.set_title("Descomposicion de la varianza (relativa)")
    ax.legend(fontsize=8)
    for i in range(len(x)):
        if ps[i] > 4:
            ax.text(i, ps[i] / 2, f"{ps[i]:.0f}%", ha="center", color="w", fontsize=8)
        if pr[i] > 4:
            ax.text(i, ps[i] + pr[i] / 2, f"{pr[i]:.0f}%", ha="center", color="w", fontsize=8)

    fig.suptitle("Las dos fuentes de varianza, SEPARADAS\n"
                 "senal = campo GP (predecible)   |   ruido = verosimilitud (irreducible)",
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_descomposicion_varianza.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_media_varianza(curvas, outdir):
    """
    LA figura que distingue las verosimilitudes: Var(X|lambda) vs lambda.
        gaussiana : plana en tau^2      (homocedastica)
        poisson   : diagonal, var = lam
        NB        : cuadratica, var = lam + lam^2/r

    curvas: {nombre_dataset: (xs, ys)} ya calculadas con
            analyze_datasets.curva_media_varianza. Los datasets que falten
            simplemente no se dibujan.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colores = dict(gaussian="seagreen", poisson="steelblue",
                   nb_s00="gray", nb_s10="darkorange", nb_s50="firebrick")

    for ax, (nombres, titulo, logy) in zip(
            axes,
            [ (["gaussian", "poisson"],
               "Gaussiana vs Poisson: la clave es si el ruido escala con la media", False),
              (["gaussian", "poisson", "nb_s10", "nb_s50"],
               "Anadiendo la NB (escala log): sobredispersion cuadratica", True) ]):
        for nm in nombres:
            if nm not in curvas:
                continue
            xs, ys = curvas[nm]
            ax.plot(xs, ys, "o-", ms=4, color=colores.get(nm), label=nm)
        lr = np.linspace(20, 90, 100)
        ax.plot(lr, lr, "k--", lw=1.2, label="var = lambda (Poisson teorico)")
        ax.axhline(50.0, color="green", ls=":", lw=1.2, label="tau^2 = 50 (gaussiano teorico)")
        ax.set_xlabel("lambda (media local)"); ax.set_ylabel("Var(X | lambda)")
        ax.set_title(titulo, fontsize=10)
        if logy:
            ax.set_yscale("log")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.suptitle("RELACION MEDIA-VARIANZA: la firma de cada verosimilitud", fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_media_varianza.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_dispersion_sparsity(lam, outdir):
    """
    Relacion ANALITICA dispersion <-> sparsity de la NB:
        P(X=0)    = (r/(r+lam))^r
        var/media = 1 + lam/r
    """
    mu = float(lam.mean())
    r_vals = np.logspace(-2, 4, 200)
    spars = np.array([np.mean((r / (r + lam))**r) for r in r_vals])
    disp = 1 + mu / r_vals

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    ax = axes[0]
    ax.loglog(r_vals, disp, "-", color="darkorange", lw=2)
    ax.axhline(1, color="steelblue", ls="--", label="Poisson (=1)")
    ax.set_xlabel("dispersion r"); ax.set_ylabel("indice de dispersion (1 + mu/r)")
    ax.set_title("Dispersion en funcion de r")
    ax.legend(); ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(100 * spars, disp, "-", color="darkorange", lw=2)
    for target, color in [(0.10, "green"), (0.50, "red")]:
        j = int(np.argmin(np.abs(spars - target)))
        ax.plot(100 * spars[j], disp[j], "o", ms=9, color=color,
                label=f"sparsity {target*100:.0f}% -> r={r_vals[j]:.3f}, disp={disp[j]:.1f}")
    ax.axhline(1, color="steelblue", ls="--", label="Poisson (sparsity ~ 0)")
    ax.set_yscale("log")
    ax.set_xlim(-2, 100)
    ax.set_xlabel("sparsity (% de ceros)"); ax.set_ylabel("indice de dispersion")
    ax.set_title("DISPERSION EN FUNCION DE LA SPARSITY (relacion analitica)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")

    fig.suptitle(f"Binomial negativa: r controla dispersion y sparsity a la vez "
                 f"(mu={mu:.0f})", fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_dispersion_vs_sparsity.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


# --------------------------------------------------------------------------- #
# CEROS: DROPOUT vs BIOLOGICO                                                  #
# --------------------------------------------------------------------------- #
def fig_ceros(res, outdir):
    r_nb = [r for r in res if r["modelo"].startswith("nb")]
    if not r_nb:
        return None
    nombres = [r["dataset"] for r in r_nb]
    x = np.arange(len(nombres))
    bio = np.array([r["sparsity_bio_pct"] for r in r_nb])
    amb = np.array([r["ceros_que_son_mascara_pct"] for r in r_nb])

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    ax.bar(x, bio, color="firebrick")
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel("% de ceros biologicos")
    ax.set_title("Ceros que produce el proceso generativo")
    for i, v in enumerate(bio):
        ax.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=9)

    ax = axes[1]
    ax.bar(x, amb, color="seagreen")
    ax.set_xticks(x); ax.set_xticklabels(nombres)
    ax.set_ylabel("% de los ceros observados que son de mascara")
    ax.set_ylim(0, 105)
    ax.set_title("AMBIGUEDAD: cuantos ceros son recuperables\n"
                 "(el resto son biologicos: no hay que tocarlos)", fontsize=10)
    for i, v in enumerate(amb):
        ax.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=9)

    fig.suptitle("Dropout vs significado biologico: cuanto mas baja la barra verde, "
                 "mas ambigua es la tarea", fontsize=12)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_tipos_de_cero.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_matriz_confusion(resultados_conf, outdir):
    """Heatmaps de las matrices de confusion, una por dataset con ceros mixtos."""
    nombres = list(resultados_conf.keys())
    n = len(nombres)
    if n == 0:
        return None
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4))
    if n == 1:
        axes = [axes]
    for ax, nombre in zip(axes, nombres):
        c = resultados_conf[nombre]
        M = np.array([[c["tp"], c["fn"]], [c["fp"], c["tn"]]])
        im = ax.imshow(M, cmap="Blues")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{M[i,j]}", ha="center", va="center",
                        fontsize=13, color="black")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["pred. dropout", "pred. biologico"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["real dropout", "real biologico"])
        ax.set_title(f"{nombre}\nprec={c['precision']:.2f}  rec={c['recall']:.2f}  "
                     f"f1={c['f1']:.2f}", fontsize=9)
    fig.suptitle("Deteccion INGENUA del origen de cada cero (basada en vecinos espaciales)\n"
                 "sirve de referencia: cualquier metodo de imputacion real deberia superarla",
                 fontsize=11)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_matriz_confusion.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


# --------------------------------------------------------------------------- #
# DOMINIOS ESPACIALES (gen_datasets3)                                          #
# --------------------------------------------------------------------------- #
def fig_dominios_mapa(coords, g_size, variantes, outdir):
    n = len(variantes)
    fig, axes = plt.subplots(2, n, figsize=(6 * n, 9))
    for j, (modo, labels) in enumerate(variantes.items()):
        ax = axes[0, j]
        im = ax.imshow(labels.reshape(g_size, g_size), origin="lower", cmap="Set2")
        ax.set_title(f"dominios -- modo={modo}", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, ticks=np.unique(labels))

        ax = axes[1, j]
        counts_d = np.bincount(labels)
        ax.bar(range(len(counts_d)), counts_d, color="steelblue")
        ax.axhline(labels.size / len(counts_d), color="red", ls="--", lw=1, label="reparto equitativo")
        ax.set_xlabel("dominio"); ax.set_ylabel("n spots")
        ax.set_title(f"tamano por dominio (total={labels.size})", fontsize=10)
        ax.legend(fontsize=8)

    fig.suptitle("Dominios espaciales discretos: libre vs atada\n"
                 "(mismo campo lambda base, mismo n_dominios)", fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_dominios_mapa.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_moran_por_dominio(resultados_por_variante, outdir, titulo=None):
    fig, ax = plt.subplots(figsize=(9, 5))
    dominios = sorted(next(iter(resultados_por_variante.values())).keys())
    n_series = len(resultados_por_variante)
    ancho = 0.8 / n_series
    x = np.arange(len(dominios))
    for i, (modo, res) in enumerate(resultados_por_variante.items()):
        ax.bar(x + i * ancho, [res.get(d, np.nan) for d in dominios],
               width=ancho, label=modo)
    ax.set_xticks(x + ancho * (n_series - 1) / 2)
    ax.set_xticklabels([f"dominio {d}" for d in dominios])
    ax.set_ylabel("Moran's I medio (sobre genes, k=8)")
    ax.set_title(titulo or "Coherencia espacial interna de cada dominio")
    ax.axhline(0, color="black", lw=0.8)
    ax.legend(); ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_moran_por_dominio.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p


def fig_correlacion_por_variante(lam_libre, lam_atada, asignacion, outdir):
    """Heatmap gen-gen lado a lado: revela si 'atada' introduce bloques de
    dominio en la matriz de correlacion (contaminacion dominio<->modulo)."""
    orden = np.argsort(asignacion)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, lam, modo in zip(axes, [lam_libre, lam_atada], ["libre", "atada"]):
        C = np.corrcoef(lam[:, orden].T)
        im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_title(f"modo={modo}", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axes, fraction=0.03, label="correlacion de Pearson")
    fig.suptitle("Correlacion gen-gen: libre vs atada", fontsize=12)
    p = os.path.join(outdir, "fig_correlacion_por_variante.png")
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p
