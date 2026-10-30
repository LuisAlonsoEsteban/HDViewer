"""Translate application settings into explicit backend arguments.

The backend has standalone defaults and never imports this module. Application
settings are applied here at object creation or analysis launch; explicit UI
choices override the defaults in these mappings.
"""
from copy import deepcopy

from . import *


def sample_kwargs(kind=None):
    options = dict(gene_cache_size=GENE_CACHE_SIZE,
                   min_gene_cells=MIN_GENE_CELLS,
                   normalization_target_sum=NORMALIZATION_TARGET_SUM,
                   marker_method=DE_METHOD)
    if kind == 'bins':
        options.update(aggregated_prefix=AGGREGATED_PREFIX,
                       aggregation_factors=AGGREGATION_FACTORS,
                       max_cache_entries=MAX_CACHE_ENTRIES,
                       aggregate_obs_prefix=AGGREGATE_OBS_PREFIX)
    elif kind == 'cells':
        options.update(lod_prefix=LOD_PREFIX, lod_base_tolerance=LOD_BASE_TOLERANCE,
                       lod_max_level=LOD_MAX_LEVEL, lod_tolerance_px=LOD_TOLERANCE_PX,
                       lod_hysteresis=LOD_HYSTERESIS)
    return options


def grid_fit_kwargs():
    return dict(max_points=GRID_FIT_MAX_POINTS, max_condition_number=GRID_MAX_CONDITION_NUMBER)


def bins_to_spots_kwargs():
    return dict(grid_fit_max_points=GRID_FIT_MAX_POINTS,
                grid_max_condition_number=GRID_MAX_CONDITION_NUMBER,
                spot_buffer_quad_segs=SPOT_BUFFER_QUAD_SEGS,
                normalization_target_sum=NORMALIZATION_TARGET_SUM)


def image_kwargs():
    return dict(tile_size=TILE_SIZE)


def clustering_kwargs():
    return dict(resolution=CLUSTER_RESOLUTION, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS,
                algorithm=CLUSTER_ALGORITHM, n_clusters=CLUSTER_N_CLUSTERS,
                eps=CLUSTER_DBSCAN_EPS, min_samples=CLUSTER_MIN_SAMPLES,
                min_cluster_size=CLUSTER_MIN_SIZE, covariance_type=CLUSTER_COVARIANCE_TYPE,
                linkage=CLUSTER_LINKAGE, random_state=CLUSTER_RANDOM_STATE)


def embedding_kwargs():
    return dict(n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)


def deconvolution_kwargs():
    return dict(batch_size=DECONV_BATCH_SIZE)


def comparison_kwargs():
    return dict(min_coverage=COMPARISON_MIN_COVERAGE, predicate=COMPARISON_PREDICATE,
                n_hvg=COMPARISON_N_HVG, n_top_svg=COMPARISON_TOP_SPATIAL_GENES,
                normalization_target_sum=NORMALIZATION_TARGET_SUM,
                chunk_size=COMPARISON_CHUNK_SIZE, k_neighbors=COMPARISON_K_NEIGHBORS,
                hull_ratio=COMPARISON_HULL_RATIO)


def _merge(base, overrides):
    result = deepcopy(base)
    for key, value in overrides.items():
        result[key] = _merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else deepcopy(value)
    return result


def segmentation_kwargs():
    from ..startup.data.image.segmentation.service import load_config
    defaults = dict(runner=SEGMENTATION_DEFAULT_RUNNER, gpu=SEGMENTATION_GPU,
                    threads_per_run=DEFAULT_THREADS, max_parallel=SEGMENTATION_MAX_PARALLEL,
                    docker=dict(binary=SEGMENTATION_DOCKER_BINARY),
                    apptainer=dict(binary=SEGMENTATION_APPTAINER_BINARY, images_dir=SEGMENTATION_IMAGES_DIR))
    return dict(config=_merge(defaults, load_config(SEGMENTATION_CONFIG_PATH)),
                config_path=SEGMENTATION_CONFIG_PATH, plugins_dir=PLUGINS_DIR,
                sdk_dir=SDK_DIR, segmentation_root=SEGMENTATION_ROOT)