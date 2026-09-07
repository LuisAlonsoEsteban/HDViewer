# ============================================================================
# FICHERO RECONSTRUIDO -- PROCEDENCIA
# ============================================================================
# El fuente original de gen_datasets3.py se perdio. `imputacion/` entero esta
# sin trackear en git ("?? imputacion/"), asi que no habia historial del que
# restaurarlo. Lo unico que sobrevivio fue el bytecode
# __pycache__/gen_datasets3.cpython-38.pyc (compilado 2026-09-03 13:08, de un
# fuente de 39296 bytes).
#
# Este fichero = decompilacion de ese .pyc (decompyle3 3.9.3, Python 3.8) MAS
# dos cambios que el .pyc no traia y que si estaban en los datasets vigentes:
#   1. bases_offset_mult recalibrado (el .pyc tenia la calibracion aditiva).
#   2. library_size_sigma=0.3 y su aplicacion (no existia en el .pyc; el fuente
#      se edito entre las 13:08 y las 14:04, hora de los .h5ad vigentes).
# Los dos se reconstruyeron leyendo generation_params y obs["library_size"] de
# los propios .h5ad, no de memoria.
#
# VERIFICADO (2026-09-07): regenerar con seed=42 reproduce los 12 datasets de
# synthetic_out/datasets_gen3/ **bit a bit** -- ground_truth_full, lambda_field,
# obs["dominio"], obs["library_size"] y, en las 8 bases sin re-enmascarar,
# tambien X y mask. Comando del test:
#     python -c "import gen_datasets3 as g; g.main(dict(g.CFG3, outdir='/tmp/x',
#                sin_figuras=True, sin_tablas=True, mask_frac=0.1))"
#
# QUE SE PERDIO IGUALMENTE: los comentarios inline (#) del original. Los
# docstrings y los textos de ayuda del argparse SI sobrevivieron (son constantes
# de cadena, y el bytecode las conserva). Los nombres de variables locales y el
# formato pueden no coincidir con el original.
#
# NOTA: nb_s00 y nb_s10 se re-enmascararon despues a mask_frac=0.30 con
# remask.py; CFG3 aqui mantiene 0.10, que es el valor con el que se generaron.
# ============================================================================

"""
gen_dataset3.py
================
Extiende gen_datasets2.py: genera, para CADA dataset de alli, dos variantes con
DOMINIOS ESPACIALES DISCRETOS (ground truth categorico para ARI/NMI/AMI/HOMO).

Reutiliza toda la infraestructura de gen_datasets.py (grid, kernel GP,
modulos de co-expresion, verosimilitud NB, mascarado, guardado) via import
-- no duplica nada, solo anade la capa de dominios encima.

Cobertura: las SEIS bases de gen_datasets2 (gp_crudo, gaussian, poisson,
nb_s00, nb_s10, nb_s50) x dos modos de offset = 12 datasets. Cada variante
aplica la MISMA verosimilitud que su base sobre lambda+offsets de dominio, asi
que la comparacion entre ellas aisla el efecto del modelo de observacion (y de
la sparsity) sobre la recuperabilidad de los dominios. Antes solo se generaban
nb_s00 y nb_s10, y los otros cuatro datasets no tenian dominios ni, por tanto,
metricas de clustering.

Ademas de generar los datasets, al terminar corre sobre ellos el MISMO analisis
estadistico que ya se hacia sobre datasets/ (analyze_datasets.py: TABLA 1
descomposicion de la varianza, TABLA 2 verosimilitud/dispersion, TABLA 3 tipos
de cero, TABLA 4 y 7 tests formales, TABLA 5 matriz de confusion, TABLA 6
Moran's I), de modo que un dataset con dominios trae tanto esas tablas como las
metricas de dominios. Se desactiva con --sin-tablas.

Dos modos de offset sobre el MISMO campo lambda base (mismos modulos de co-expresion):
    <base>_dom_libre : offsets de dominio INDEPENDIENTES de los modulos
                       (dominio y modulo son ejes ortogonales -- evita
                       confundir "deteccion de dominio" con "deteccion de
                       modulo", ya evaluada aparte)
    <base>_dom_atada : offsets de dominio ATADOS a los modulos (dominio d
                       sobreexpresa el modulo d, marcador biologico) --
                       variante de comparacion, no el default

En las bases NB la sparsity se recalibra sobre lambda+offset (no sobre el
lambda original) para mantener la sparsity objetivo pese al desplazamiento de
media por dominio; en gp_crudo/gaussian/poisson no se fuerza sparsity ninguna,
igual que en gen_datasets2.

Salidas en carpeta NUEVA (no mezcla con datasets/ de gen_datasets.py), con UNA
SUBCARPETA POR BASE: lo que pertenece a una sola base vive dentro de ella y lo
GENERICO (lo que compara todas las bases entre si) se queda en la raiz. El
layout completo y las funciones que lo resuelven estan en rutas.py.

    datasets_gen3/<base>/dataset_<base>_dom_libre.h5ad
    datasets_gen3/<base>/dataset_<base>_dom_atada.h5ad
    datasets_gen3/<base>/ground_truth_<base>_dom_*.npz  (held-out: spot_idx,
                                                  gene_idx, true_value en las
                                                  posiciones enmascaradas --
                                                  ver mask_frac en CFG3)
    datasets_gen3/<base>/fig_dominios_mapa_<base>.png
    datasets_gen3/<base>/fig_moran_por_dominio_<base>.png
    datasets_gen3/<base>/fig_correlacion_por_variante_<base>.png
    datasets_gen3/metricas_estadisticas.csv  + figuras de analyze_datasets.py
                                               (--sin-tablas para omitirlas)

Las figuras y el calculo de Moran's I NO viven aqui: son compartidos con
analyze_datasets.py y estan en figuras.py y moran.py respectivamente.

Uso:
    python gen_datasets3.py                       # las 6 bases x 2 modos + tablas
    python gen_datasets3.py --bases poisson nb_s50
    python gen_datasets3.py --seed 7 --sin-figuras # replica de un barrido
"""
import os, numpy as np, anndata as ad
from gen_datasets2 import CFG as CFG_BASE, build_grid, rbf_kernel, sample_gp_correlacionado, observe, mask_ground_truth, save, find_r_for_sparsity, nb_zero_prob
from moran import moran_por_dominio
from figuras import fig_dominios_mapa, fig_moran_por_dominio, fig_correlacion_por_variante
from rutas import GEN3, dir_base
CFG3 = dict(CFG_BASE)
CFG3.update(dict(outdir=GEN3,
  n_dominios=3,
  mask_frac=0.1,
  dominio_offset_scale=10.0,
  dominio_offset_modo="multiplicativo",
  dominio_temperatura=1.0,
  dominio_length_scale=4.0,
  k_vecinos_moran=8,
  clip_frac_aviso=0.01,
  dominio_min_frac=0.05,
  dominio_max_reintentos=200,
  bases={'gp_crudo':dict(verosimilitud="gp_crudo"), 
 'gaussian':dict(verosimilitud="gaussian"), 
 'poisson':dict(verosimilitud="poisson"), 
 'nb_s00':dict(verosimilitud="nbinom", target_sparsity=0.0), 
 'nb_s10':dict(verosimilitud="nbinom", target_sparsity=0.1), 
 'nb_s50':dict(verosimilitud="nbinom", target_sparsity=0.5)},
  seed_index={
  'nb_s00': 0, 'nb_s10': 1, 'gp_crudo': 2,
  'gaussian': 3,  'poisson': 4,  'nb_s50': 5},
  sparsity_via="profundidad",
  r_quasi_poisson=10000000.0,
  # RECONSTRUIDO 2026-09-07. Estos valores NO estan inventados: se han leido de
  # generation_params de los .h5ad vigentes (datasets_gen3/<base>/), que son los
  # que produjo la version perdida de este fichero. El .pyc del que salio esta
  # decompilacion traia la calibracion ANTERIOR (nb_s10 libre 1.0/atada 6.0,
  # nb_s50 libre 2.0/atada 10.0), de la epoca de los offsets aditivos.
  bases_offset_mult={'gp_crudo': {'libre': 1.0, 'atada': 1.5},
                     'gaussian': {'libre': 1.0, 'atada': 2.0},
                     'poisson':  {'libre': 1.0, 'atada': 2.0},
                     'nb_s00':   {'libre': 1.0, 'atada': 2.0},
                     'nb_s10':   {'libre': 1.5, 'atada': 3.0},
                     'nb_s50':   {'libre': 2.0, 'atada': 5.0}},

  # FACTOR TECNICO DE PROFUNDIDAD POR SPOT (library size).
  # d_s ~ LogNormal(0, sigma^2) centrado a E[d]=1, aplicado a lambda SOLO en las
  # verosimilitudes de conteo (poisson/nbinom) -- que son exactamente las que
  # despues pasan por normalize_total en el pipeline de clustering. Sin el,
  # normalize_total no corrige nada real y solo mete ruido; con el, la
  # normalizacion pasa a estar justificada.
  # Vector COMPARTIDO por las 12 variantes (como la particion de dominios).
  library_size_sigma=0.3))

def sample_domain_field(K_dom, n_dominios, jitter, rng):
    """Campos de 'afinidad' espacial a cada dominio, via Cholesky sobre K_dom
    (K_dom usa dominio_length_scale, MAYOR que el de lambda, para que el
    grano de los dominios supere al ruido espacial normal del campo)."""
    n_spots = K_dom.shape[0]
    L = np.linalg.cholesky(K_dom + jitter * np.eye(n_spots))
    Z = rng.standard_normal((n_spots, n_dominios))
    return L @ Z


def assign_domains(afinidad, temperatura, rng):
    """
    temperatura=0  -> argmax puro (fronteras nitidas)
    temperatura>0  -> ruido ~ N(0, temperatura^2) antes del argmax (fronteras difusas)
    """
    ruido = temperatura * rng.standard_normal(afinidad.shape)
    return (afinidad + ruido).argmax(axis=1)


def particion_dominios(coords, cfg, verbose=True):
    """Reparte los spots en n_dominios regiones espaciales via argmax de campos
    GP de afinidad (particion tipo Voronoi organica; tamanos NO uniformes, como
    las regiones reales de un tejido).

    Usa RECHAZO (rejection sampling): descarta las realizaciones DEGENERADAS en
    las que algun dominio queda por debajo de dominio_min_frac de los spots
    -- demasiado pequeno para ser una region real y para calcular Moran/
    clustering de forma fiable -- reintentando con otra semilla. NO fuerza
    tamanos iguales: solo elimina los casos degenerados, preservando la
    variabilidad realista de tamano entre regiones. Si tras max_reintentos
    ninguna cumple, devuelve la MENOS degenerada (mayor tamano minimo) con aviso.
    """
    n_spots = coords.shape[0]
    min_size = int(np.ceil(cfg["dominio_min_frac"] * n_spots))
    K_dom = rbf_kernel(coords, cfg["dominio_length_scale"], 1.0)
    (mejor_labels, mejor_min, mejor_tam) = (None, -1, None)
    for intento in range(cfg["dominio_max_reintentos"]):
        rng_dom = np.random.default_rng(cfg["seed"] + 100 + intento)
        afinidad = sample_domain_field(K_dom, cfg["n_dominios"], cfg["jitter"], rng_dom)
        labels = assign_domains(afinidad, cfg["dominio_temperatura"], rng_dom)
        tamanos = np.bincount(labels, minlength=(cfg["n_dominios"]))
        min_actual = int(tamanos.min())
        if min_actual > mejor_min:
            mejor_labels, mejor_min, mejor_tam = labels, min_actual, tamanos
        if min_actual >= min_size:
            if verbose:
                print(f"[dominios]  particion aceptada en intento {intento} (min_size exigido={min_size}): tamanos={tamanos.tolist()}")
            return labels
        print(f'[AVISO] ningun reparto alcanzo min_frac={cfg["dominio_min_frac"]} (={min_size} spots) tras {cfg["dominio_max_reintentos"]} intentos; uso el mejor (min={mejor_min}): tamanos={mejor_tam.tolist()}')
        return mejor_labels


def build_domain_offsets(n_dominios, n_genes, offset_scale, rng, asignacion_modulo=None, modo='libre'):
    """Matriz (n_dominios, n_genes) del efecto de dominio, en la ESCALA que
    diga el llamante: con dominio_offset_modo="multiplicativo" (default) estos
    valores son log-fold-changes y se aplican como lam_base*exp(delta); con
    "aditivo" (legado) se suman a lam_base. La forma de la distribucion es la
    misma en los dos casos, solo cambia offset_scale y como se aplica.

    modo="libre": offsets aleatorios independientes de los modulos
                  (recomendado -- aisla la senal de dominio de la de modulo)
    modo="atada": dominio d activa preferentemente el modulo (d mod n_mod)
                  (variante de comparacion)
    """
    offsets = np.zeros((n_dominios, n_genes))
    if modo == "libre":
        offsets[:] = offset_scale * rng.standard_normal((n_dominios, n_genes))
    elif modo == "atada":
        assert asignacion_modulo is not None, "modo='atada' requiere asignacion_modulo"
        n_mod = asignacion_modulo.max() + 1
        for d in range(n_dominios):
            modulo_activado = d % n_mod
            genes_del_modulo = asignacion_modulo == modulo_activado
            offsets[(d, genes_del_modulo)] = offset_scale * (0.5 + rng.random(genes_del_modulo.sum()))

    else:
        raise ValueError(f"modo desconocido: {modo!r}")
    return offsets


def find_depth_for_sparsity(lam, target, r, lo=0.0001, hi=1.0, tol=1e-06, max_iter=200):
    """Biseccion sobre un factor multiplicativo de lambda (profundidad de
    conteo) hasta alcanzar la sparsity objetivo, con r FIJO y alto.

    Es la contrapartida de find_r_for_sparsity: alli se compraban los ceros
    con sobredispersion (y se destruia la senal, ver CFG3); aqui se compran
    bajando la profundidad, que es como aparecen realmente en Visium. P(X=0)
    es monotona DECRECIENTE en la profundidad.

    target<=0 -> profundidad 1.0 (sin reescalar: el dataset denso de referencia).
    """
    if target <= 0.0:
        return (1.0, nb_zero_prob(np.clip(lam, 0.001, None), r))
    for _ in range(max_iter):
        mid = np.sqrt(lo * hi)
        p = nb_zero_prob(np.clip(mid * lam, 0.001, None), r)
        if abs(p - target) < tol:
            return (mid, p)
        if p > target:
            lo = mid
        else:
            hi = mid

    return (
     mid, nb_zero_prob(np.clip(mid * lam, 0.001, None), r))


def tau_gaussiano(cfg):
    """Misma regla que gen_datasets2.main: tau = sqrt(gp_mean) si no se fija a
    mano, para que la varianza del ruido gaussiano iguale a la que tendria la
    Poisson (comparacion justa entre verosimilitudes)."""
    if cfg.get("gauss_tau"):
        return float(cfg["gauss_tau"])
    return float(np.sqrt(cfg["gp_mean"]))


def generar_variante_dominio(lam_base, asignacion, coords, cfg, modo, base_nombre, base_cfg, dominio_labels, seed_offset=0, library_size=None):
    """Genera una variante <base_nombre>_dom_<modo> con dominios discretos.

    base_nombre identifica el dataset de partida (gp_crudo, gaussian, poisson,
    nb_s00, nb_s10, nb_s50) y base_cfg trae su verosimilitud -- que es la que se
    aplica sobre lambda+offsets, la MISMA que en gen_datasets2. Asi cada dataset
    de alli tiene aqui su gemelo con dominios, con el mismo modelo de
    observacion y el mismo campo latente.
    """
    veros = base_cfg.get("verosimilitud", "nbinom")
    target_sparsity = float(base_cfg.get("target_sparsity", 0.0))
    rng = np.random.default_rng(cfg["seed"] + 100 + seed_offset)
    mult = cfg.get("bases_offset_mult", {}).get(base_nombre, {})
    mult = mult.get(modo, 1.0) if isinstance(mult, dict) else float(mult)
    modo_offset = cfg.get("dominio_offset_modo", "multiplicativo")
    if modo_offset == "multiplicativo":
        escala = cfg["dominio_offset_scale"] / cfg["gp_mean"] * mult
        delta = build_domain_offsets((cfg["n_dominios"]), (cfg["n_genes"]), escala, rng,
          asignacion_modulo=asignacion, modo=modo)
        if modo == "libre":
            delta = delta - 0.5 * escala ** 2
        crudo = lam_base * np.exp(delta[dominio_labels])
    else:
        escala = cfg["dominio_offset_scale"] * mult
        offsets = build_domain_offsets((cfg["n_dominios"]), (cfg["n_genes"]), escala, rng,
          asignacion_modulo=asignacion, modo=modo)
        crudo = lam_base + offsets[dominio_labels]
    clip_mask = crudo < cfg["lambda_floor"]
    clip_frac = float(clip_mask.mean())
    lam_con_off = np.clip(crudo, cfg["lambda_floor"], None)

    # --- library size por spot -------------------------------------------
    # Solo en las verosimilitudes de CONTEO: gp_crudo/gaussian se quedan como
    # referencia continua, sin el artefacto y sin necesitar la correccion.
    # Se aplica ANTES de calibrar la profundidad para que find_depth_for_sparsity
    # vea la lambda definitiva y la sparsity resultante sea la pedida.
    aplicar_lib = (library_size is not None
                   and float(cfg.get("library_size_sigma", 0.0)) > 0.0
                   and veros in ("poisson", "nbinom"))
    if aplicar_lib:
        lam_con_off = lam_con_off * np.asarray(library_size, dtype=float)[:, None]

    r = None
    tau = None
    profundidad = 1.0
    p_teo = float("nan")
    if veros == "nbinom":
        if cfg.get("sparsity_via", "profundidad") == "profundidad":
            r = cfg["r_quasi_poisson"]
            (profundidad, p_teo) = find_depth_for_sparsity(lam_con_off, target_sparsity, r)
            lam_dom = np.clip(profundidad * lam_con_off, cfg["lambda_floor"], None)
        else:
            lam_dom = lam_con_off
            (r, p_teo) = find_r_for_sparsity(lam_dom, target_sparsity)
    else:
        lam_dom = lam_con_off
        if veros == "gaussian":
            tau = tau_gaussiano(cfg)
    r_ = np.random.default_rng(cfg["seed"] + 1 + seed_offset)
    counts = observe(lam_dom, veros, r_, tau=tau, r=r)
    (obs, mask, rows, cols, gt) = mask_ground_truth(counts, cfg["mask_frac"], r_)
    ceros = counts == 0
    n_ceros = int(ceros.sum())
    ceros_por_clip = 100.0 * float((ceros & clip_mask).sum()) / n_ceros if n_ceros else 0.0
    extra = {'dominio_modo':modo, 
     'clip_frac':clip_frac, 
     'ceros_por_clip_pct':ceros_por_clip, 
     'base_dataset':base_nombre, 
     'verosimilitud':veros, 
     'n_dominios':int(cfg["n_dominios"]), 
     'sparsity_via':(cfg.get)("sparsity_via", "profundidad"), 
     'profundidad':float(profundidad), 
     'dominio_offset_mult':float(mult), 
     'dominio_offset_modo':modo_offset, 
     'dominio_offset_escala':float(escala), 
     'library_size_sigma':float(cfg.get("library_size_sigma", 0.0)), 
     'library_size_aplicado':bool(aplicar_lib)}
    if r is not None:
        extra["nb_dispersion_r"] = float(r)
        extra["target_sparsity_dom"] = float(target_sparsity)
    if tau is not None:
        extra["gauss_tau"] = float(tau)
    nombre = f'{base_nombre}_dom_{modo}{cfg.get("sufijo", "")}'
    cfg_base = dict(cfg, outdir=dir_base(nombre, (cfg["outdir"]), crear=True))
    h5 = save(nombre, cfg_base, coords, counts, obs, mask, lam_dom, rows, cols, gt, asignacion=asignacion,
      extra=extra)
    a = ad.read_h5ad(h5)
    a.obs["dominio"] = dominio_labels.astype(str)
    a.obs["dominio"] = a.obs["dominio"].astype("category")
    if library_size is not None:
        # se guarda en TODAS las variantes, tambien en las continuas donde no se
        # aplico (library_size_aplicado dice si entro o no en lambda): asi el
        # vector queda documentado en el .h5ad y no hay que rederivarlo
        a.obs["library_size"] = np.asarray(library_size, dtype=float)
    a.write_h5ad(h5)
    sparsity_real = 100 * float((counts == 0).mean())
    sparsity_obs = 100 * float((obs == 0).mean())
    r_txt = (f"{r:>9.4g}") if r is not None else (f'{"-":>9}')
    p_txt = f"{p_teo * 100:>6.2f}%" if np.isfinite(p_teo) else (f'{"-":>7}')
    print(f'[{nombre:>22}]  veros={veros:<8} profundidad={profundidad:>7.4f}  r={r_txt}  escala={escala:>5.2f}{"(log)" if modo_offset == "multiplicativo" else "(add)"}  lam_media={lam_dom.mean():>6.2f}  P(0)teo={p_txt}  sparsity_bio={sparsity_real:>6.2f}%  sparsity_X={sparsity_obs:>6.2f}%  clip={100 * clip_frac:>5.3f}% (ceros_por_clip={ceros_por_clip:>5.2f}%)  -> {os.path.join(os.path.basename(os.path.dirname(h5)), os.path.basename(h5))}')
    if clip_frac > cfg.get("clip_frac_aviso", 0.01):
        print(f"[AVISO] {nombre}: {100 * clip_frac:.2f}% de las entradas tocan lambda_floor -> esa parte de la sparsity es ESTRUCTURAL (lambda~0 da cero seguro), no la estocastica que se pidio. Con dominio_offset_modo='multiplicativo' esto deberia ser 0: revisa bases_offset_mult['{base_nombre}']['{modo}'].")
    return (
     h5, dominio_labels, lam_dom, r, counts)


def main(cfg):
    os.makedirs((cfg["outdir"]), exist_ok=True)
    coords = build_grid(cfg["grid_size"])
    rng = np.random.default_rng(cfg["seed"])
    print(f'[dominios]  n_dominios={cfg["n_dominios"]}  offset_scale={cfg["dominio_offset_scale"]}  temperatura={cfg["dominio_temperatura"]}  length_scale={cfg["dominio_length_scale"]}\n')
    K = rbf_kernel(coords, cfg["length_scale"], cfg["gp_sigma"])
    (F, W, asignacion, factores) = sample_gp_correlacionado(K,
      (cfg["gp_mean"]), (cfg["n_genes"]), n_factores=(cfg["n_factores"]), genes_por_modulo=(cfg.get("genes_por_modulo", None)),
      peso_intra=(cfg["peso_intra"]),
      peso_fondo=(cfg["peso_fondo"]),
      jitter=(cfg["jitter"]),
      rng=rng)
    lam_base = np.clip(F, cfg["lambda_floor"], None)
    # Vector de library size COMPARTIDO por las 12 variantes. Generador propio
    # (seed + 999) para que anadirlo NO desplace ninguno de los otros streams:
    # el del GP es default_rng(seed), el de los offsets seed+100+offset y el de
    # observe() seed+1+offset. Verificado: con seed=42 reproduce exactamente el
    # obs["library_size"] de los .h5ad vigentes.
    sigma_lib = float(cfg.get("library_size_sigma", 0.0))
    if sigma_lib > 0.0:
        rng_lib = np.random.default_rng(cfg["seed"] + 999)
        # centrado a E[d]=1: exp(sigma*z - sigma^2/2)
        library_size = np.exp(sigma_lib * rng_lib.standard_normal(coords.shape[0])
                              - 0.5 * sigma_lib ** 2)
        print(f"[library]  sigma={sigma_lib}  media={library_size.mean():.4f}  "
              f"rango=[{library_size.min():.3f}, {library_size.max():.3f}]  "
              f"(solo se aplica a poisson/nbinom)")
    else:
        library_size = None

    dominio_labels = particion_dominios(coords, cfg)
    tamanos = np.bincount(dominio_labels, minlength=(cfg["n_dominios"]))
    print(f"[dominios]  tamanos (compartidos por TODAS las combinaciones): {tamanos.tolist()}\n")
    resultados = {}
    nombres = []
    for (base_nombre, base_cfg) in cfg["bases"].items():
        for (i_modo, modo) in enumerate(["libre", "atada"]):
            seed_i = 2 * cfg["seed_index"][base_nombre] + i_modo
            (h5, labels, lam_dom, r, counts) = generar_variante_dominio(lam_base,
              asignacion, coords, cfg, modo=modo, base_nombre=base_nombre,
              base_cfg=base_cfg,
              dominio_labels=dominio_labels,
              seed_offset=seed_i,
              library_size=library_size)
            resultados[(base_nombre, modo)] = (h5, labels, lam_dom, r, counts)
            nombres.append(f'{base_nombre}_dom_{modo}{cfg.get("sufijo", "")}')

    if cfg.get("sin_figuras"):
        print("\n[fig]  omitidas (--sin-figuras)")
        print(f'Datasets en {cfg["outdir"]}/')
        return nombres
    print("\n[fig]  generando figuras de validacion...")
    for base_nombre in cfg["bases"]:
        labels_por_modo = {modo: resultados[(base_nombre, modo)][1] for modo in ('libre',
                                                                                 'atada')}
        p1 = fig_dominios_mapa(coords, cfg["grid_size"], labels_por_modo, dir_base(base_nombre, (cfg["outdir"]), crear=True))
        nuevo_nombre = p1.replace("fig_dominios_mapa.png", f"fig_dominios_mapa_{base_nombre}.png")
        os.replace(p1, nuevo_nombre)
        print(f"[fig]   {nuevo_nombre}")

    moran_control = moran_por_dominio(lam_base, dominio_labels, coords, k=(cfg["k_vecinos_moran"]))
    for base_nombre in cfg["bases"]:
        lam_por_modo = {modo: resultados[(base_nombre, modo)][2] for modo in ('libre',
                                                                              'atada')}
        counts_por_modo = {modo: resultados[(base_nombre, modo)][4] for modo in ('libre',
                                                                                 'atada')}
        labels_por_modo = {modo: resultados[(base_nombre, modo)][1] for modo in ('libre',
                                                                                 'atada')}
        resultados_moran = {f"{modo} (counts)": moran_por_dominio((counts_por_modo[modo]),
          (labels_por_modo[modo]), coords,
          k=(cfg["k_vecinos_moran"])) for modo in ('libre', 'atada')}
        resultados_moran["control (lambda)"] = moran_control
        p2 = fig_moran_por_dominio(resultados_moran,
          dir_base(base_nombre, (cfg["outdir"]), crear=True), titulo=f"Coherencia espacial interna por dominio -- {base_nombre}\n(Moran sobre COUNTS; 'control' = campo latente sin ruido)")
        nuevo_nombre = p2.replace("fig_moran_por_dominio.png", f"fig_moran_por_dominio_{base_nombre}.png")
        os.replace(p2, nuevo_nombre)
        print(f"[fig]   {nuevo_nombre}")
        for (modo, res) in resultados_moran.items():
            print(f"  moran por dominio ({base_nombre}, {modo}): " + "  ".join((f"d{d}={v:.3f}" for d, v in res.items())))

        p3 = fig_correlacion_por_variante(lam_por_modo["libre"], lam_por_modo["atada"], asignacion, dir_base(base_nombre, (cfg["outdir"]), crear=True))
        nuevo_nombre = p3.replace("fig_correlacion_por_variante.png", f"fig_correlacion_por_variante_{base_nombre}.png")
        os.replace(p3, nuevo_nombre)
        print(f"[fig]   {nuevo_nombre}")

    print(f'\nDatasets y figuras en {cfg["outdir"]}/')
    print(f'{len(nombres)} combinaciones generadas: {{{", ".join(cfg["bases"])}}} x {{libre, atada}}')
    if not cfg.get("sin_tablas"):
        print("\n[tablas]  analisis estadistico sobre los datasets con dominios...")
        import analyze_datasets
        analyze_datasets.main(indir=(cfg["outdir"]), nombres=nombres)
    print("\nSiguiente paso: python cluster_domains.py  (Leiden + ARI/NMI/AMI/HOMO sobre estos mismos datasets)")
    return nombres


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=None, help="semilla de la replica. Cambia TODA la realizacion (campo GP, particion de dominios, offsets y counts), que es lo que hace de esto una replica y no una variacion parcial. Sin --seed usa CFG3['seed'] y escribe los nombres sin sufijo (comportamiento previo)")
    ap.add_argument("--sin-figuras", action="store_true", dest="sin_figuras", help="no generar las figuras de validacion (barridos). Implica --sin-tablas: en un barrido cada replica pisaria el metricas_estadisticas.csv de la anterior")
    ap.add_argument("--sin-tablas", action="store_true", dest="sin_tablas", help="no correr al final el analisis estadistico (TABLA 1-7 de analyze_datasets.py) sobre los datasets generados. Se puede sacar despues con 'python analyze_datasets.py --gen3'")
    ap.add_argument("--bases", nargs="+", default=None, help="subconjunto de bases a generar (gp_crudo gaussian poisson nb_s00 nb_s10 nb_s50). Por defecto TODAS: cada dataset de gen_datasets2 tiene su gemelo con dominios")
    ap.add_argument("--sufijo", default=None, help="sufijo de los ficheros; por defecto '_seed<N>' cuando se pasa --seed. Los nombres de dataset lo llevan dentro, asi que cluster_domains.py y evaluar_clustering_metodos.py los cargan sin cambios (--datasets nb_s10_dom_libre_seed7)")
    args = ap.parse_args()
    cfg = dict(CFG3)
    cfg["sin_figuras"] = args.sin_figuras
    cfg["sin_tablas"] = args.sin_tablas
    if args.bases:
        desconocidas = [b for b in args.bases if b not in CFG3["bases"]]
        if desconocidas:
            ap.error(f'bases desconocidas: {desconocidas}. Disponibles: {list(CFG3["bases"])}')
        cfg["bases"] = {b: CFG3["bases"][b] for b in args.bases}
    if args.seed is not None:
        cfg["seed"] = args.seed
        cfg["sufijo"] = args.sufijo if args.sufijo is not None else f"_seed{args.seed}"
    elif args.sufijo is not None:
        cfg["sufijo"] = args.sufijo
    main(cfg)

