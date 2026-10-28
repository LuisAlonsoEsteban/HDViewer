"""Translate application settings into explicit backend arguments.

The backend has standalone defaults and never imports this module. Application
settings are applied here at object creation or analysis launch; explicit UI
choices override the defaults in these mappings.
"""
from copy import deepcopy

from . import settings


def sample_kwargs(kind=None):
    options = dict(gene_cache_size=settings.GENE_CACHE_SIZE, cluster_col=settings.CLUSTER_COL,
                   min_gene_cells=settings.MIN_GENE_CELLS,
                   normalization_target_sum=settings.NORMALIZATION_TARGET_SUM,
                   marker_method=settings.DE_METHOD)
    if kind == 'bins':
        options.update(aggregated_prefix=settings.AGGREGATED_PREFIX,
                       aggregation_factors=settings.AGGREGATION_FACTORS,
                       max_cache_entries=settings.MAX_CACHE_ENTRIES,
                       cluster_prefix=settings.CLUSTER_PREFIX,
                       aggregate_obs_prefix=settings.AGGREGATE_OBS_PREFIX)
    elif kind == 'cells':
        options.update(lod_prefix=settings.LOD_PREFIX, lod_base_tolerance=settings.LOD_BASE_TOLERANCE,
                       lod_max_level=settings.LOD_MAX_LEVEL, lod_tolerance_px=settings.LOD_TOLERANCE_PX,
                       lod_hysteresis=settings.LOD_HYSTERESIS)
    return options


def grid_fit_kwargs():
    return dict(max_points=settings.GRID_FIT_MAX_POINTS, max_condition_number=settings.GRID_MAX_CONDITION_NUMBER)


def bins_to_spots_kwargs():
    return dict(grid_fit_max_points=settings.GRID_FIT_MAX_POINTS,
                grid_max_condition_number=settings.GRID_MAX_CONDITION_NUMBER,
                spot_buffer_quad_segs=settings.SPOT_BUFFER_QUAD_SEGS,
                normalization_target_sum=settings.NORMALIZATION_TARGET_SUM)


def image_kwargs():
    return dict(tile_size=settings.TILE_SIZE)


def clustering_kwargs():
    return dict(resolution=settings.CLUSTER_RESOLUTION, n_pcs=settings.N_PCS, n_neighbors=settings.N_NEIGHBORS,
                algorithm=settings.CLUSTER_ALGORITHM, n_clusters=settings.CLUSTER_N_CLUSTERS,
                eps=settings.CLUSTER_DBSCAN_EPS, min_samples=settings.CLUSTER_MIN_SAMPLES,
                min_cluster_size=settings.CLUSTER_MIN_SIZE, covariance_type=settings.CLUSTER_COVARIANCE_TYPE,
                linkage=settings.CLUSTER_LINKAGE, random_state=settings.CLUSTER_RANDOM_STATE)


def embedding_kwargs():
    return dict(n_pcs=settings.N_PCS, n_neighbors=settings.N_NEIGHBORS)


def deconvolution_kwargs():
    return dict(batch_size=settings.DECONV_BATCH_SIZE)


def comparison_kwargs():
    return dict(min_coverage=settings.COMPARISON_MIN_COVERAGE, predicate=settings.COMPARISON_PREDICATE,
                n_hvg=settings.COMPARISON_N_HVG, n_top_svg=settings.COMPARISON_TOP_SPATIAL_GENES,
                normalization_target_sum=settings.NORMALIZATION_TARGET_SUM,
                chunk_size=settings.COMPARISON_CHUNK_SIZE, k_neighbors=settings.COMPARISON_K_NEIGHBORS,
                hull_ratio=settings.COMPARISON_HULL_RATIO)


def _merge(base, overrides):
    result = deepcopy(base)
    for key, value in overrides.items():
        result[key] = _merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else deepcopy(value)
    return result


def segmentation_kwargs():
    from src.segmentation.service import load_config
    defaults = dict(runner=settings.SEGMENTATION_DEFAULT_RUNNER, gpu=settings.SEGMENTATION_GPU,
                    threads_per_run=settings.DEFAULT_THREADS, max_parallel=settings.SEGMENTATION_MAX_PARALLEL,
                    docker=dict(binary=settings.SEGMENTATION_DOCKER_BINARY),
                    apptainer=dict(binary=settings.SEGMENTATION_APPTAINER_BINARY, images_dir=settings.SEGMENTATION_IMAGES_DIR))
    return dict(config=_merge(defaults, load_config(settings.SEGMENTATION_CONFIG_PATH)),
                config_path=settings.SEGMENTATION_CONFIG_PATH, plugins_dir=settings.PLUGINS_DIR,
                sdk_dir=settings.SDK_DIR, segmentation_root=settings.SEGMENTATION_ROOT)
