"""Translate application settings into explicit backend arguments.

The backend has standalone defaults and never imports this module. Application
settings are applied here at object creation or analysis launch; explicit UI
choices override the defaults in these mappings.
"""
from copy import deepcopy

from .. import settings as _settings


def sample_kwargs(kind=None):
    options = dict(gene_cache_size=_settings.GENE_CACHE_SIZE,
                   min_gene_cells=_settings.MIN_GENE_CELLS,
                   normalization_target_sum=_settings.NORMALIZATION_TARGET_SUM,
                   marker_method=_settings.DE_METHOD)
    if kind == 'bins':
        options.update(aggregated_prefix=_settings.AGGREGATED_PREFIX,
                       aggregation_factors=_settings.AGGREGATION_FACTORS,
                       max_cache_entries=_settings.MAX_CACHE_ENTRIES,
                       aggregate_obs_prefix=_settings.AGGREGATE_OBS_PREFIX)
    elif kind == 'cells':
        options.update(lod_prefix=_settings.LOD_PREFIX, lod_base_tolerance=_settings.LOD_BASE_TOLERANCE,
                       lod_max_level=_settings.LOD_MAX_LEVEL, lod_tolerance_px=_settings.LOD_TOLERANCE_PX,
                       lod_hysteresis=_settings.LOD_HYSTERESIS)
    return options


def grid_fit_kwargs():
    return dict(max_points=_settings.GRID_FIT_MAX_POINTS, max_condition_number=_settings.GRID_MAX_CONDITION_NUMBER)


def bins_to_spots_kwargs():
    return dict(grid_fit_max_points=_settings.GRID_FIT_MAX_POINTS,
                grid_max_condition_number=_settings.GRID_MAX_CONDITION_NUMBER,
                spot_buffer_quad_segs=_settings.SPOT_BUFFER_QUAD_SEGS,
                normalization_target_sum=_settings.NORMALIZATION_TARGET_SUM)


def image_kwargs():
    return dict(tile_size=_settings.TILE_SIZE)


def clustering_kwargs():
    return dict(resolution=_settings.CLUSTER_RESOLUTION, n_pcs=_settings.N_PCS, n_neighbors=_settings.N_NEIGHBORS,
                algorithm=_settings.CLUSTER_ALGORITHM, n_clusters=_settings.CLUSTER_N_CLUSTERS,
                eps=_settings.CLUSTER_DBSCAN_EPS, min_samples=_settings.CLUSTER_MIN_SAMPLES,
                min_cluster_size=_settings.CLUSTER_MIN_SIZE, covariance_type=_settings.CLUSTER_COVARIANCE_TYPE,
                linkage=_settings.CLUSTER_LINKAGE, random_state=_settings.CLUSTER_RANDOM_STATE)


def embedding_kwargs():
    return dict(n_pcs=_settings.N_PCS, n_neighbors=_settings.N_NEIGHBORS)


def deconvolution_kwargs():
    return dict(batch_size=_settings.DECONV_BATCH_SIZE)


def comparison_kwargs():
    return dict(min_coverage=_settings.COMPARISON_MIN_COVERAGE, predicate=_settings.COMPARISON_PREDICATE,
                n_hvg=_settings.COMPARISON_N_HVG, n_top_svg=_settings.COMPARISON_TOP_SPATIAL_GENES,
                normalization_target_sum=_settings.NORMALIZATION_TARGET_SUM,
                chunk_size=_settings.COMPARISON_CHUNK_SIZE, k_neighbors=_settings.COMPARISON_K_NEIGHBORS,
                hull_ratio=_settings.COMPARISON_HULL_RATIO)


def _merge(base, overrides):
    result = deepcopy(base)
    for key, value in overrides.items():
        result[key] = _merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else deepcopy(value)
    return result


def segmentation_kwargs():
    from ..startup.data.image.segmentation.service import load_config
    defaults = dict(runner=_settings.SEGMENTATION_DEFAULT_RUNNER, gpu=_settings.SEGMENTATION_GPU,
                    threads_per_run=_settings.DEFAULT_THREADS, max_parallel=_settings.SEGMENTATION_MAX_PARALLEL,
                    timeout_minutes=_settings.SEGMENTATION_TIMEOUT_MINUTES, memory=_settings.SEGMENTATION_MEMORY,
                    pids_limit=_settings.SEGMENTATION_PIDS_LIMIT,
                    docker=dict(binary=_settings.SEGMENTATION_DOCKER_BINARY),
                    apptainer=dict(binary=_settings.SEGMENTATION_APPTAINER_BINARY, images_dir=_settings.SEGMENTATION_IMAGES_DIR))
    return dict(config=_merge(defaults, load_config(_settings.SEGMENTATION_CONFIG_PATH)),
                config_path=_settings.SEGMENTATION_CONFIG_PATH, plugins_dir=_settings.PLUGINS_DIR,
                sdk_dir=_settings.SDK_DIR, segmentation_root=_settings.SEGMENTATION_ROOT)


def deconvolution_service_kwargs():
    from ..board.operations.layer.deconvolution.service import load_config
    defaults = dict(runner=_settings.DECONVOLUTION_DEFAULT_RUNNER, gpu=_settings.DECONVOLUTION_GPU,
                    threads_per_run=_settings.DEFAULT_THREADS, max_parallel=_settings.DECONVOLUTION_MAX_PARALLEL,
                    timeout_minutes=_settings.DECONVOLUTION_TIMEOUT_MINUTES, memory=_settings.DECONVOLUTION_MEMORY,
                    pids_limit=_settings.DECONVOLUTION_PIDS_LIMIT,
                    work_dir=str(_settings.DECONVOLUTION_WORK_DIR))
    return dict(config=_merge(defaults, load_config(_settings.DECONVOLUTION_CONFIG_PATH)),
                root=_settings.DECONVOLUTION_ROOT, plugins_dir=_settings.DECONVOLUTION_ROOT / 'plugins',
                sdk_dir=_settings.DECONVOLUTION_ROOT / 'sdk')
