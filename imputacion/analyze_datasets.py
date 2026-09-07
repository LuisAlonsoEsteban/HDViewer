r"""
analyze_datasets.py
=====================
SCRIPT 2 de 2: analisis estadistico de los datasets de gen_datasets.py.

Responde a las tres preguntas del tutor:

  (1) "no solapes el ruido gaussiano con el de Poisson"
      -> TABLA 1: descomposicion de la varianza. Separa numericamente las dos
         fuentes usando la ley de la varianza total:
              Var(X) = Var(E[X|lam])  +  E[Var(X|lam)]
                       \_ del GP ___/    \_ de la verosimilitud _/
                       = la SENAL         = el RUIDO

  (2) "que dataset es mejor"
      -> TABLA 2 + fig_media_varianza: compara las verosimilitudes por su
         relacion media-varianza, que es lo que las distingue de verdad:
              gaussiana : Var(X|lam) = tau^2        CONSTANTE (homocedastica)
              poisson   : Var(X|lam) = lam          crece con la media
              NB        : Var(X|lam) = lam + lam^2/r  crece cuadraticamente

  (3) "distinguir si los 0s son dropout o tienen significado biologico"
      -> TABLA 3 + fig_tipos_de_cero: cuantifica la ambiguedad para cada nivel
         de sparsity (0%, 10%, 50%).

  (4) "estadisticos importantes: p-value, matriz de confusion"
      -> TABLA 4: test formal de sobredispersion (p-value) -- es significativamente
         distinto de Poisson, o no.
      -> TABLA 5 + fig_matriz_confusion: matriz de confusion de un detector
         INGENUO de dropout (basado en vecinos espaciales), con precision/recall.
         Es el mismo principio que el baseline kNN de evaluate_imputation.py,
         aplicado aqui a la pregunta de clasificar el ORIGEN de cada cero.

Este script se queda con los CALCULOS (descomposicion, dispersion, tests) y el
main que imprime las tablas. Lo que estaba aqui y era compartido vive ahora en
dos modulos aparte:
    moran.py    -> autocorrelacion espacial (Moran's I y su test)
    figuras.py  -> todas las figuras (las de aqui y las de gen_datasets3.py)

Uso:
    python analyze_datasets.py
"""

import re
import numpy as np
import anndata as ad
from scipy import stats
import os
import csv

from rutas import DATASETS as INDIR, GEN3, ruta_dataset, datasets_en
from moran import build_knn_weights, test_moran_permutacion
from figuras import (
    fig_descomposicion, fig_media_varianza, fig_dispersion_sparsity,
    fig_ceros, fig_matriz_confusion,
)

ORDEN = ["gp_crudo", "gaussian", "poisson", "nb_s00", "nb_s10", "nb_s50"]

# los datasets con dominios de gen_datasets3.py: las mismas bases, cada una en
# sus dos modos de offset. Se analizan con ESTE mismo codigo (la capa de
# dominios no cambia ni la verosimilitud ni el campo lambda guardado), por eso
# main() solo necesita saber de que carpeta leer.
ORDEN_GEN3 = [f"{b}_dom_{m}" for b in ORDEN for m in ("libre", "atada")]


def verosimilitud_de(params, nombre):
    """Familia de observacion del dataset: 'gp_crudo', 'gaussian', 'poisson' o
    'nbinom'. Es lo que decide que test aplica cada tabla.

    gen_datasets3.py la guarda EXPLICITA en generation_params['verosimilitud'],
    porque alli 'modelo' es el nombre completo del dataset
    (gaussian_dom_libre, nb_s10_dom_atada) y ya no identifica la familia. En
    los datasets de gen_datasets2 el nombre SI es la familia -> fallback al
    nombre, que es lo que este script uso siempre.
    """
    v = params.get("verosimilitud")
    if v:
        return str(v)
    return str(params.get("modelo", nombre))


# --------------------------------------------------------------------------- #
# (1) DESCOMPOSICION DE LA VARIANZA                                            #
# --------------------------------------------------------------------------- #
def descomponer_varianza(counts, lam, modelo, tau=None, r=None):
    """
    Ley de la varianza total:
        Var(X) = Var(E[X|lam])  +  E[Var(X|lam)]
                 \___senal___/     \____ruido____/

    Var(E[X|lam]) = Var(lam) = sigma^2   -> la SENAL espacial, comun a todos
    E[Var(X|lam)] depende de la verosimilitud:
        gp_crudo : 0                     (no hay ruido de observacion)
        gaussian : tau^2
        poisson  : E[lam] = mu
        nbinom   : E[lam + lam^2/r] = mu + E[lam^2]/r
    """
    var_senal = float(lam.var())
    mu = float(lam.mean())

    if modelo == "gp_crudo":
        var_ruido = 0.0
    elif modelo == "gaussian":
        var_ruido = float(tau)**2
    elif modelo == "poisson":
        var_ruido = mu
    elif modelo.startswith("nb"):
        var_ruido = mu + float((lam**2).mean()) / float(r)
    else:
        var_ruido = np.nan

    var_teo = var_senal + var_ruido
    total = var_teo if var_teo > 0 else np.nan
    return dict(
        var_empirica=float(counts.var()),
        var_teorica=var_teo,
        var_senal=var_senal,
        var_ruido=var_ruido,
        pct_senal=100.0 * var_senal / total if total else np.nan,
        pct_ruido=100.0 * var_ruido / total if total else np.nan,
        snr=var_senal / var_ruido if var_ruido > 0 else np.inf,
    )


# --------------------------------------------------------------------------- #
# (2) DISPERSION: MARGINAL vs CONDICIONAL                                      #
# --------------------------------------------------------------------------- #
def dispersion_condicional(counts, lam, n_bins=20):
    """
    var/media agrupando por bins de lambda parecido: compara 'manzanas con
    manzanas'. Es la dispersion DE LA VEROSIMILITUD, sin contaminar con la
    variacion del campo. Poisson -> ~1.  NB -> ~1 + lambda/r.
    """
    lam_f, c_f = lam.ravel(), counts.ravel()
    bins = np.quantile(lam_f, np.linspace(0, 1, n_bins + 1))
    idx = np.digitize(lam_f, bins)
    ratios = []
    for b in range(1, len(bins)):
        m = idx == b
        if m.sum() > 10 and c_f[m].mean() > 0:
            ratios.append(c_f[m].var() / c_f[m].mean())
    return float(np.mean(ratios)) if ratios else np.nan


def curva_media_varianza(counts, lam, n_bins=18):
    """
    Var(X|lam) empirica en funcion de lam. ESTA es la firma de cada
    verosimilitud: gaussiana plana, Poisson diagonal, NB cuadratica.
    """
    lam_f, c_f = lam.ravel(), counts.ravel()
    bins = np.quantile(lam_f, np.linspace(0, 1, n_bins + 1))
    idx = np.digitize(lam_f, bins)
    xs, ys = [], []
    for b in range(1, len(bins)):
        m = idx == b
        if m.sum() > 10:
            xs.append(lam_f[m].mean())
            ys.append(c_f[m].var())
    return np.array(xs), np.array(ys)


# --------------------------------------------------------------------------- #
# ANALISIS DE UN DATASET                                                       #
# --------------------------------------------------------------------------- #
def analizar(nombre, indir=INDIR):
    a = ad.read_h5ad(ruta_dataset(nombre, indir))
    p = dict(a.uns["generation_params"])
    counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
    lam = np.asarray(a.layers["lambda_field"], dtype=float)
    mask = np.asarray(a.layers["mask"]).astype(bool)
    obs = np.asarray(a.X, dtype=float)

    modelo = verosimilitud_de(p, nombre)
    tau = p.get("gauss_tau")
    r = p.get("nb_dispersion_r")

    vd = descomponer_varianza(counts, lam, modelo, tau=tau, r=r)

    n_bio = int((counts == 0).sum())
    n_obs_zero = int((obs == 0).sum())
    n_mask_zero = int((mask & (counts > 0)).sum())
    media = float(counts.mean())
    es_conteo = modelo not in ("gp_crudo", "gaussian")

    return dict(
        dataset=nombre, modelo=modelo,
        # en gen3 el nombre del dataset ya no dice la verosimilitud: se guarda
        # aparte para que el CSV siga siendo legible por si solo
        dominio=str(p.get("dominio_modo", "-")),
        r=(float(r) if r is not None else np.nan),
        media=media,
        # (1) descomposicion
        var_emp=vd["var_empirica"], var_teo=vd["var_teorica"],
        var_senal=vd["var_senal"], var_ruido=vd["var_ruido"],
        pct_senal=vd["pct_senal"], pct_ruido=vd["pct_ruido"], snr=vd["snr"],
        # (2) dispersion + validez
        disp_marginal=(vd["var_empirica"] / media if media > 0 else np.nan),
        disp_condicional=(dispersion_condicional(counts, lam) if es_conteo else np.nan),
        disp_cond_teo=((1 + float(lam.mean()) / float(r)) if r is not None
                       else (1.0 if modelo == "poisson" else np.nan)),
        pct_negativos=100.0 * float((counts < 0).mean()),
        es_entero=bool(np.all(counts == np.round(counts))),
        # (3) ceros
        sparsity_bio_pct=100.0 * n_bio / counts.size,
        sparsity_obs_pct=100.0 * n_obs_zero / counts.size,
        n_ceros_bio=n_bio, n_mascara=int(mask.sum()),
        ceros_que_son_mascara_pct=(100.0 * n_mask_zero / n_obs_zero
                                   if n_obs_zero else np.nan),
        gt_que_son_cero_pct=(100.0 * float((counts[mask] == 0).mean())
                             if mask.sum() else np.nan),
    )


# --------------------------------------------------------------------------- #
# (4a) TEST FORMAL DE SOBREDISPERSION (con p-value)                            #
# --------------------------------------------------------------------------- #
def test_sobredispersion(counts, lam):
    """
    Test de Cameron & Trivedi (1990) para equidispersion de Poisson.

    H0: los datos son Poisson (Var(X|lam) = lam)
    H1: los datos estan sobredispersos (Var(X|lam) > lam)

    Estadistico (usa lambda VERDADERO como valor ajustado, porque aqui se
    conoce exactamente -- en datos reales se usaria la media ajustada de un
    GLM Poisson):
        z = sum( (x_i - lam_i)^2 - x_i )  /  sqrt( 2 * sum(lam_i^2) )
    Bajo H0, z ~ N(0,1) asintoticamente. Test de una cola (la sobredispersion
    va en una unica direccion).

    Devuelve (z, p_value). p_value pequeno -> se rechaza Poisson pura.
    """
    x = counts.ravel(); l = lam.ravel()
    numer = np.sum((x - l)**2 - x)
    denom = np.sqrt(2.0 * np.sum(l**2))
    z = numer / denom
    p = 1.0 - stats.norm.cdf(z)     # una cola: H1 = sobredispersion (z grande)
    return float(z), float(p)


def test_chi2_bondad_ajuste(counts, lam, modelo, r=None, tau=None):
    """
    Test chi-cuadrado de bondad de ajuste: dentro del bin de lambda mas
    cercano a la media global, compara el histograma EMPIRICO de counts
    contra las frecuencias TEORICAS del modelo (Poisson, NB o Normal).

    H0: los datos siguen la distribucion teorica del modelo.
    p-value pequeno -> se rechaza esa distribucion como generadora de los datos.

    Los bins se agrupan de forma ADAPTATIVA (acumulando categorias hasta que
    la frecuencia esperada >= 5, el requisito estandar del test) formando una
    PARTICION COMPLETA del rango observado -- asi la suma observada y la
    esperada coinciden exactamente, sin truncar cola a mano.
    """
    mu_ref = float(lam.mean())
    lam_f, x_f = lam.ravel(), counts.ravel()
    cerca = np.abs(lam_f - mu_ref) < max(2.0, 0.05 * mu_ref)
    x_sub = x_f[cerca]
    n = len(x_sub)
    if n < 30:
        return np.nan, np.nan, n

    if modelo in ("gp_crudo", "gaussian"):
        if n > 5000:
            x_sub = np.random.default_rng(0).choice(x_sub, 5000, replace=False)
        stat, p = stats.shapiro(x_sub)
        return float(stat), float(p), n

    if modelo == "poisson" or (modelo.startswith("nb") and r is not None and r > 1e5):
        dist = stats.poisson(mu_ref)
    elif modelo.startswith("nb"):
        dist = stats.nbinom(r, r / (r + mu_ref))
    else:
        return np.nan, np.nan, n

    k_max = int(x_sub.max())
    ks = np.arange(0, k_max + 1)
    pmf = dist.pmf(ks)
    obs_counts = np.array([np.sum(x_sub == k) for k in ks], dtype=float)
    tail_prob = max(0.0, 1.0 - pmf.sum())   # masa teorica por encima de k_max

    grupos_obs, grupos_exp = [], []
    cur_o, cur_e = 0.0, 0.0
    for o_k, p_k in zip(obs_counts, pmf):
        cur_o += o_k; cur_e += p_k * n
        if cur_e >= 5:
            grupos_obs.append(cur_o); grupos_exp.append(cur_e)
            cur_o, cur_e = 0.0, 0.0
    # cola final (lo que no llego a 5 + la masa teorica > k_max): se funde
    # con el ultimo grupo para no perder ninguna observacion (particion completa)
    cur_e += tail_prob * n
    if grupos_exp:
        grupos_obs[-1] += cur_o
        grupos_exp[-1] += cur_e
    else:
        grupos_obs, grupos_exp = [cur_o], [cur_e]

    if len(grupos_exp) < 2:
        return np.nan, np.nan, n

    obs_arr = np.array(grupos_obs, dtype=float)
    exp_arr = np.array(grupos_exp, dtype=float)
    exp_arr *= obs_arr.sum() / exp_arr.sum()   # corrige redondeo de punto flotante
    chi2, p = stats.chisquare(obs_arr, exp_arr)
    return float(chi2), float(p), n


# --------------------------------------------------------------------------- #
# (4d) TEST DE LEVENE: homogeneidad de varianza entre bins de lambda           #
# --------------------------------------------------------------------------- #
def test_levene_homocedasticidad(counts, lam, n_bins=6):
    """
    H0: la varianza de X es la MISMA en todos los bins de lambda
        (homocedasticidad -- la firma de la verosimilitud GAUSSIANA).
    H1: la varianza CAMBIA segun lambda
        (heterocedasticidad -- la firma de Poisson/NB, donde Var crece con la media).

    Formaliza con un p-value lo que se ve en fig_media_varianza.
    p < 0.05 -> se rechaza homocedasticidad (varianza NO constante).
    """
    lam_f, x_f = lam.ravel(), counts.ravel()
    bins = np.quantile(lam_f, np.linspace(0, 1, n_bins + 1))
    idx = np.digitize(lam_f, bins)
    grupos = [x_f[idx == b] for b in range(1, n_bins + 1) if (idx == b).sum() > 10]
    if len(grupos) < 2:
        return np.nan, np.nan
    stat, p = stats.levene(*grupos, center="median")
    return float(stat), float(p)


# --------------------------------------------------------------------------- #
# (4e) VERIFICACION: eliminar el artefacto de agrupar por "lambda parecido"    #
# --------------------------------------------------------------------------- #
# Las TABLA 4 originales agrupan por bins de lambda APROXIMADO (cuantiles).
# Eso mezcla observaciones con lambda ligeramente distintos dentro del mismo
# bin, lo que puede inflar artificialmente tanto el chi2 (poisson rechazo H0
# con p=0.021) como el Levene (gaussian salio "heterocedastico" cuando no
# deberia). Aqui se COMPRUEBA esa hipotesis usando el lambda EXACTO de cada
# observacion individual -- posible solo porque el dataset es sintetico y se
# conoce el lambda real sin error. Si el artefacto es la explicacion correcta,
# estos tests deberian dejar de rechazar donde antes rechazaban.

def test_pit_uniformidad(counts, lam, modelo, r=None, tau=None, seed=0):
    """
    PIT (Probability Integral Transform) + Kolmogorov-Smirnov, usando el
    lambda EXACTO de CADA observacion (sin ningun binning por "parecido").

    Si el modelo es correcto, u_i = F(x_i; lambda_i) ~ Uniforme(0,1),
    para CUALQUIER valor de lambda_i -- no hace falta agrupar nada.
    Para distribuciones DISCRETAS (Poisson, NB) se usa el PIT ALEATORIZADO
    (Dunn & Smyth, 1996) para no romper la uniformidad por los saltos de
    la CDF discreta:
        u_i ~ Uniform( F(x_i - 1; lambda_i), F(x_i; lambda_i) )
    """
    rng = np.random.default_rng(seed)
    x, l = counts.ravel(), lam.ravel()

    if modelo == "gp_crudo":
        return np.nan, np.nan          # sin ruido de observacion: PIT no aplica
    if modelo == "gaussian":
        u = stats.norm.cdf(x, loc=l, scale=tau)
    elif modelo == "poisson":
        u = rng.uniform(stats.poisson.cdf(x - 1, l), stats.poisson.cdf(x, l))
    elif modelo.startswith("nb"):
        p_fail = r / (r + l)
        u = rng.uniform(stats.nbinom.cdf(x - 1, r, p_fail),
                        stats.nbinom.cdf(x, r, p_fail))
    else:
        return np.nan, np.nan

    stat, p = stats.kstest(u, "uniform")
    return float(stat), float(p)


def test_levene_residuos_exactos(counts, lam, n_bins=6):
    """
    Repite Levene pero sobre los RESIDUOS (counts - lambda_EXACTO de cada
    punto) en vez de sobre los counts crudos, agrupados por los MISMOS bins
    de lambda que antes.

    Si la heterocedasticidad "detectada" antes en 'gaussian' era un artefacto
    de mezclar lambdas ligeramente distintos dentro de cada bin (en vez de una
    propiedad real del ruido), al restar el lambda EXACTO ese artefacto
    desaparece: el test deberia dejar de ser significativo para 'gaussian'.
    Para poisson/NB, en cambio, el residuo AUN depende de lambda (su varianza
    es lambda o lambda+lambda^2/r), asi que seguira siendo significativo --
    eso es correcto, no es el artefacto.
    """
    resid = counts - lam
    lam_f, r_f = lam.ravel(), resid.ravel()
    bins = np.quantile(lam_f, np.linspace(0, 1, n_bins + 1))
    idx = np.digitize(lam_f, bins)
    grupos = [r_f[idx == b] for b in range(1, n_bins + 1) if (idx == b).sum() > 10]
    if len(grupos) < 2:
        return np.nan, np.nan
    stat, p = stats.levene(*grupos, center="median")
    return float(stat), float(p)


# --------------------------------------------------------------------------- #
# (4b) MATRIZ DE CONFUSION: detector INGENUO del origen de cada cero           #
# --------------------------------------------------------------------------- #
def detector_ingenuo_dropout(obs, mask, counts, coords, k=8):
    """
    Para cada entrada observada como CERO, predice si es "dropout" (recuperable)
    o "biologico" (de verdad no se expresa), usando SOLO informacion disponible
    para un metodo real (nunca la mascara ni el ground truth):

        prediccion = "dropout"   si MAS DE LA MITAD de los k vecinos espaciales
                                  NO enmascarados expresan ese gen (voto mayoritario)
        prediccion = "biologico" en caso contrario

    (Un criterio mas laxo como "algun vecino lo expresa" predice "dropout" casi
    siempre y no discrimina nada -- por eso se usa voto mayoritario.)

    Es el mismo principio que el baseline kNN de evaluate_imputation.py,
    aplicado aqui como CLASIFICADOR en vez de como regresor.

    Verdad (solo conocida porque el dataset es sintetico):
        "dropout"   = mask=True  y el valor real (antes de enmascarar) > 0
        "biologico" = valor real == 0  (fue enmascarada o no, da igual: seguiria
                      siendo 0 de todas formas)

    Devuelve matriz de confusion + metricas.
    """
    n_spots = coords.shape[0]
    d2 = ((coords[:, None, :] - coords[None, :, :])**2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    orden = np.argsort(d2, axis=1)

    ceros_i, ceros_g = np.where(obs == 0)
    if len(ceros_i) == 0:
        return None

    y_true, y_pred = [], []
    for i, g in zip(ceros_i, ceros_g):
        es_dropout_real = bool(mask[i, g] and counts[i, g] > 0)
        vecinos = orden[i][:k]
        # voto mayoritario: dropout si MAS DE LA MITAD de los k vecinos
        # (contando los enmascarados como "no observado", se ignoran) expresan
        # el gen. Un criterio de "algun vecino" es demasiado laxo -> predice
        # dropout casi siempre y no discrimina nada.
        obs_vecinos = obs[vecinos, g]
        no_enmascarados = ~mask[vecinos, g]
        votos_nonzero = np.sum((obs_vecinos > 0) & no_enmascarados)
        votos_total = np.sum(no_enmascarados)
        pred_dropout = (votos_total > 0) and (votos_nonzero / votos_total > 0.5)
        y_true.append(es_dropout_real)
        y_pred.append(pred_dropout)

    y_true = np.array(y_true); y_pred = np.array(y_pred)
    tp = int(np.sum(y_true & y_pred))
    fp = int(np.sum(~y_true & y_pred))
    fn = int(np.sum(y_true & ~y_pred))
    tn = int(np.sum(~y_true & ~y_pred))

    precision = tp / (tp + fp) if (tp + fp) else np.nan
    recall = tp / (tp + fn) if (tp + fn) else np.nan
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and (precision + recall) > 0 else np.nan)
    accuracy = (tp + tn) / len(y_true)

    return dict(tp=tp, fp=fp, fn=fn, tn=tn, n=len(y_true),
                precision=precision, recall=recall, f1=f1, accuracy=accuracy)


# --------------------------------------------------------------------------- #
# MAIN                                                                         #
# --------------------------------------------------------------------------- #
def main(indir=INDIR, nombres=None, con_semillas=False):
    """Imprime las tablas 1-7 sobre los datasets de una carpeta.

    indir=DATASETS (default) -> los de gen_datasets2.py
    indir=GEN3               -> los de gen_datasets3.py, que son los mismos
                                modelos de observacion con dominios espaciales
                                encima; el analisis no cambia porque la capa
                                lambda_field ya incluye el offset de dominio.
    nombres: orden explicito de datasets; por defecto el ORDEN conocido de esa
             carpeta mas cualquier otro dataset_*.h5ad que haya, para no
             dejarse ninguno fuera.
    con_semillas: incluir las replicas <nombre>_seed<N> del barrido. Por
             defecto NO: son otra pregunta (variabilidad entre realizaciones,
             ver agregar_semillas.py) y triplicarian el largo de las tablas.
             Mismo criterio que cluster_domains.datasets_disponibles().
    """
    if nombres is None:
        conocidos = ORDEN_GEN3 if os.path.abspath(indir) == os.path.abspath(GEN3) else ORDEN
        # datasets_en mira tambien dentro de las subcarpetas por base de
        # datasets_gen3/ (ver rutas.py); en datasets/, que es plano, da igual
        extra = datasets_en(indir)
        if not con_semillas:
            extra = [n for n in extra if not re.search(r"_seed\d+$", n)]
        nombres = [n for n in conocidos if n in extra] + \
                  [n for n in extra if n not in conocidos]
    disp = [n for n in nombres if os.path.exists(ruta_dataset(n, indir))]
    if not disp:
        raise SystemExit(f"No hay datasets en {indir}/. "
                         f"Corre antes: python gen_datasets2.py (o gen_datasets3.py)")

    # los nombres de gen3 (nb_s10_dom_atada) no caben en la columna de 11 que
    # servia para gen2 (nb_s10): el ancho se calcula, no se fija
    w = max(11, max(len(n) for n in disp) + 2)
    ancho = 100 + (w - 11)

    res = [analizar(n, indir) for n in disp]

    # ---------------- TABLA 1 ----------------
    print("=" * ancho)
    print("TABLA 1 -- DESCOMPOSICION DE LA VARIANZA")
    print("           Var(X) = Var(lam) + E[Var(X|lam)]  =  SENAL + RUIDO")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'media':>8}{'var_emp':>11}{'var_teo':>11}"
          f"{'var_SENAL':>11}{'var_RUIDO':>11}{'%senal':>9}{'%ruido':>9}{'SNR':>9}")
    print("-" * ancho)
    for r in res:
        snr = "inf" if np.isinf(r["snr"]) else f"{r['snr']:.3f}"
        print(f"{r['dataset']:<{w}}{r['media']:>8.2f}{r['var_emp']:>11.2f}"
              f"{r['var_teo']:>11.2f}{r['var_senal']:>11.2f}{r['var_ruido']:>11.2f}"
              f"{r['pct_senal']:>9.1f}{r['pct_ruido']:>9.1f}{snr:>9}")
    print("=" * ancho)
    print("  var_emp vs var_teo: si cuadran, el modelo generativo es correcto.")
    print("  SENAL = varianza del campo GP -> es la MISMA en todos (mismo lambda).")
    print("  RUIDO = varianza que anade la verosimilitud -> es lo UNICO que cambia.")
    print("  Esto es 'no solapar': cada dataset tiene UNA sola fuente de ruido.")

    # ---------------- TABLA 2 ----------------
    print("\n" + "=" * ancho)
    print("TABLA 2 -- VEROSIMILITUD: dispersion y validez del modelo")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'r':>11}{'disp_MARG':>11}{'disp_COND':>11}{'(teo)':>10}"
          f"{'negativos%':>12}{'entero?':>10}")
    print("-" * ancho)
    for r in res:
        rs = f"{r['r']:.4g}" if not np.isnan(r["r"]) else "-"
        dc = f"{r['disp_condicional']:.3f}" if not np.isnan(r["disp_condicional"]) else "-"
        dt = f"{r['disp_cond_teo']:.3f}" if not np.isnan(r["disp_cond_teo"]) else "-"
        print(f"{r['dataset']:<{w}}{rs:>11}{r['disp_marginal']:>11.3f}{dc:>11}{dt:>10}"
              f"{r['pct_negativos']:>12.3f}{str(r['es_entero']):>10}")
    print("=" * ancho)
    print("  disp_MARG: var/media de TODA la matriz. MEZCLA senal y ruido -> para la")
    print("             Poisson NO vale 1, vale 1 + sigma^2/mu. NO es un error.")
    print("  disp_COND: var/media agrupando por lambda parecido. ESTA es la de la")
    print("             verosimilitud: ~1 para Poisson (correcto), 1+lam/r para NB.")
    print("  negativos% / entero?: la gaussiana puede dar valores negativos y no")
    print("             enteros -> imposibles en datos de conteo. Con mu=50 no ocurre")
    print("             (Poisson(50) ~ N(50,50) por el TCL), pero si con mu baja.")

    # ---------------- TABLA 3 ----------------
    print("\n" + "=" * ancho)
    print("TABLA 3 -- LOS DOS TIPOS DE CERO (dropout vs significado biologico)")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'spars_bio%':>12}{'spars_obs%':>12}{'n_ceros_bio':>13}"
          f"{'n_mascara':>11}{'ceros_mask%':>13}{'gt_cero%':>10}")
    print("-" * ancho)
    for r in res:
        cm = (f"{r['ceros_que_son_mascara_pct']:.1f}"
              if not np.isnan(r["ceros_que_son_mascara_pct"]) else "-")
        print(f"{r['dataset']:<{w}}{r['sparsity_bio_pct']:>12.3f}"
              f"{r['sparsity_obs_pct']:>12.3f}{r['n_ceros_bio']:>13}"
              f"{r['n_mascara']:>11}{cm:>13}{r['gt_que_son_cero_pct']:>10.2f}")
    print("=" * ancho)
    print("  spars_bio% : ceros del PROCESO -> significado BIOLOGICO")
    print("  ceros_mask%: de los ceros que VE el metodo, cuantos son de mascara")
    print("               (recuperables). El resto son biologicos: NO tocarlos.")
    print("               CUANTO MAS BAJO, MAS AMBIGUA la tarea de imputacion.")
    print("  gt_cero%   : % del ground truth que valia 0 (el metodo debe acertar")
    print("               esos tambien: predecir 0, no 'rellenar')")

    # ---------------- TABLA 4: tests formales (p-values) ----------------
    print("\n" + "=" * ancho)
    print("TABLA 4 -- TESTS ESTADISTICOS FORMALES")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'test':<28}{'estadistico':>13}{'p_value':>12}  {'conclusion':<28}")
    print("-" * ancho)
    for n in disp:
        a = ad.read_h5ad(ruta_dataset(n, indir))
        p = dict(a.uns["generation_params"])
        modelo = verosimilitud_de(p, n)
        counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
        lam = np.asarray(a.layers["lambda_field"], dtype=float)
        tau = p.get("gauss_tau"); r = p.get("nb_dispersion_r")

        # (a) sobredispersion (Cameron-Trivedi), solo tiene sentido en conteos
        if modelo not in ("gp_crudo", "gaussian"):
            z, pv = test_sobredispersion(counts, lam)
            concl = "sobredisperso" if pv < 0.05 else "compatible con Poisson"
            print(f"{n:<{w}}{'sobredispersion (C-T)':<28}{z:>13.3f}{pv:>12.4g}  {concl:<28}")

        # (b) bondad de ajuste
        stat, pv, nn = test_chi2_bondad_ajuste(counts, lam, modelo, r=r, tau=tau)
        if not np.isnan(pv):
            nombre_test = "Shapiro-Wilk (normalidad)" if modelo in ("gp_crudo", "gaussian") \
                          else "chi2 bondad de ajuste"
            concl = "rechaza H0" if pv < 0.05 else "no rechaza H0"
            print(f"{n:<{w}}{nombre_test:<28}{stat:>13.3f}{pv:>12.4g}  {concl:<28}  (n={nn})")

        # (c) Levene: homogeneidad de varianza entre bins de lambda
        stat, pv = test_levene_homocedasticidad(counts, lam)
        if not np.isnan(pv):
            concl = "heterocedastico" if pv < 0.05 else "homocedastico"
            print(f"{n:<{w}}{'Levene (homocedasticidad)':<28}{stat:>13.3f}{pv:>12.4g}  {concl:<28}")
    print("=" * ancho)
    print("  sobredispersion: H0 = Poisson puro. p<0.05 -> sobredisperso (esperado en NB).")
    print("  bondad de ajuste: H0 = los datos siguen esa distribucion teorica.")
    print("  Levene: H0 = varianza constante entre bins de lambda (homocedastico,")
    print("          esperado en 'gaussian'). p<0.05 -> varianza depende de lambda")
    print("          (heterocedastico, esperado en poisson/NB).")

    # ---------------- TABLA 7: VERIFICACION del artefacto de binning ----------------
    print("\n" + "=" * ancho)
    print("TABLA 7 -- VERIFICACION: eliminando el artefacto de agrupar por 'lambda parecido'")
    print("           (usa el lambda EXACTO de cada punto, sin ningun binning)")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'test':<32}{'estadistico':>13}{'p_value':>12}  {'conclusion':<24}")
    print("-" * ancho)
    for n in disp:
        a = ad.read_h5ad(ruta_dataset(n, indir))
        p = dict(a.uns["generation_params"])
        modelo = verosimilitud_de(p, n)
        counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
        lam = np.asarray(a.layers["lambda_field"], dtype=float)
        tau = p.get("gauss_tau"); r = p.get("nb_dispersion_r")

        stat, pv = test_pit_uniformidad(counts, lam, modelo, r=r, tau=tau)
        if not np.isnan(pv):
            concl = "rechaza H0" if pv < 0.05 else "no rechaza H0"
            print(f"{n:<{w}}{'PIT+KS (sin binning)':<32}{stat:>13.4f}{pv:>12.4g}  {concl:<24}")

        if modelo in ("gaussian", "poisson") or modelo.startswith("nb"):
            stat, pv = test_levene_residuos_exactos(counts, lam)
            if not np.isnan(pv):
                concl = "heterocedastico" if pv < 0.05 else "homocedastico"
                print(f"{n:<{w}}{'Levene sobre RESIDUOS exactos':<32}{stat:>13.3f}{pv:>12.4g}  {concl:<24}")
    print("=" * ancho)
    print("  PIT+KS: usa el lambda EXACTO por punto (Dunn & Smyth, 1996). Si sigue")
    print("          rechazando aqui, el problema es real, NO un artefacto de binning.")
    print("  Levene/residuos: compara con la fila equivalente de la TABLA 4. Si")
    print("          'gaussian' deja de rechazar aqui, CONFIRMA que la heterocedasticidad")
    print("          de la Tabla 4 era artefacto de mezclar lambdas dentro de cada bin.")

    # ---------------- TABLA 5: matriz de confusion ----------------
    print("\n" + "=" * ancho)
    print("TABLA 5 -- MATRIZ DE CONFUSION: deteccion INGENUA del origen de cada cero")
    print("           (basada solo en vecinos espaciales, sin usar la mascara)")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'TP':>6}{'FP':>6}{'FN':>6}{'TN':>6}{'precision':>11}"
          f"{'recall':>9}{'f1':>7}{'accuracy':>10}")
    print("-" * ancho)
    resultados_conf = {}
    coords = np.asarray(ad.read_h5ad(ruta_dataset(disp[0], indir)).obsm["spatial"])
    for n in disp:
        a = ad.read_h5ad(ruta_dataset(n, indir))
        obs = np.asarray(a.X, dtype=float)
        mask = np.asarray(a.layers["mask"]).astype(bool)
        counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
        c = detector_ingenuo_dropout(obs, mask, counts, coords, k=8)
        if c is None:
            continue
        resultados_conf[n] = c
        print(f"{n:<{w}}{c['tp']:>6}{c['fp']:>6}{c['fn']:>6}{c['tn']:>6}"
              f"{c['precision']:>11.3f}{c['recall']:>9.3f}{c['f1']:>7.3f}{c['accuracy']:>10.3f}")
    print("=" * ancho)
    print("  TP=dropout real y detectado | FP=biologico pero detectado como dropout")
    print("  FN=dropout real no detectado | TN=biologico correctamente identificado")
    print("  Detector: predice 'dropout' si MAS DE LA MITAD de los k=8 vecinos")
    print("            espaciales (no enmascarados) expresan ese gen (voto mayoritario).")
    print("            Es un baseline (misma idea que el kNN de evaluate_imputation.py);")
    print("            cualquier metodo real de imputacion deberia superarlo.")

    # ---------------- TABLA 6: Moran's I ----------------
    print("\n" + "=" * ancho)
    print("TABLA 6 -- MORAN'S I: autocorrelacion espacial (test de permutacion)")
    print("=" * ancho)
    print(f"{'dataset':<{w}}{'Moran I (media genes)':>23}{'p_value':>12}  {'conclusion':<28}")
    print("-" * ancho)
    W = build_knn_weights(coords, k=8)
    lam_ref = np.asarray(ad.read_h5ad(ruta_dataset(disp[0], indir))
                         .layers["lambda_field"], dtype=float)
    I_lam, p_lam = test_moran_permutacion(lam_ref, W, n_perm=200)
    concl = "autocorrelacion SIGNIFICATIVA" if p_lam < 0.05 else "no significativa"
    print(f"{'lambda (senal pura)':<{w}}{I_lam:>23.4f}{p_lam:>12.4g}  {concl:<28}")
    for n in disp:
        a = ad.read_h5ad(ruta_dataset(n, indir))
        counts = np.asarray(a.layers["ground_truth_full"], dtype=float)
        I_obs, pv = test_moran_permutacion(counts, W, n_perm=200)
        concl = "autocorrelacion SIGNIFICATIVA" if pv < 0.05 else "no significativa"
        print(f"{n:<{w}}{I_obs:>23.4f}{pv:>12.4g}  {concl:<28}")
    print("=" * ancho)
    print("  Moran's I in [-1, 1]: >0 = spots vecinos se parecen (lo que predice el")
    print("  kernel RBF). p-value via permutacion (H0 = sin autocorrelacion espacial).")
    print("  Comparar I entre datasets: cuanto mas ruido anade la verosimilitud, mas")
    print("  BAJA I -- cuantifica cuanta estructura espacial SOBREVIVE al ruido.")

    # ---------------- CSV + figuras ----------------
    # el CSV y las figuras de abajo comparan TODOS los datasets entre si: son
    # genericas, van en la raiz de la carpeta y no en la subcarpeta de una base
    path = os.path.join(indir, "metricas_estadisticas.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader(); w.writerows(res)
    print(f"\n[save] {path}")

    # la curva teorica dispersion-vs-sparsity solo necesita UN campo lambda de
    # referencia (es el mismo en todos los datasets de una carpeta salvo por el
    # offset de dominio): se coge el de poisson si esta, y si no el primero
    # disponible -- en gen3 el dataset se llama poisson_dom_libre, no poisson
    ref = next((n for n in disp if n.startswith("poisson")), disp[0])
    a = ad.read_h5ad(ruta_dataset(ref, indir))
    lam = np.asarray(a.layers["lambda_field"], dtype=float)

    # las curvas media-varianza se calculan aqui (leer h5ad + curva_media_varianza
    # son CALCULO, no dibujo) y se pasan ya hechas a la figura
    curvas = {}
    for nm in disp:
        a_nm = ad.read_h5ad(ruta_dataset(nm, indir))
        curvas[nm] = curva_media_varianza(
            np.asarray(a_nm.layers["ground_truth_full"], dtype=float),
            np.asarray(a_nm.layers["lambda_field"], dtype=float))

    for p in [fig_descomposicion(res, indir),
              fig_media_varianza(curvas, indir),
              fig_dispersion_sparsity(lam, indir),
              fig_ceros(res, indir),
              fig_matriz_confusion(resultados_conf, indir)]:
        if p:
            print(f"[fig]  {p}")
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--gen3", action="store_true",
                   help="analizar los datasets CON DOMINIOS de gen_datasets3.py "
                        "(datasets_gen3/) en vez de los de gen_datasets2.py")
    g.add_argument("--indir", default=None,
                   help="carpeta de datasets a analizar (alternativa a --gen3)")
    ap.add_argument("--datasets", nargs="+", default=None,
                    help="analizar solo estos datasets, en este orden")
    ap.add_argument("--con-semillas", action="store_true", dest="con_semillas",
                    help="incluir tambien las replicas _seedN del barrido "
                         "(por defecto se omiten)")
    args = ap.parse_args()

    carpeta = GEN3 if args.gen3 else (args.indir or INDIR)
    main(indir=carpeta, nombres=args.datasets, con_semillas=args.con_semillas)
