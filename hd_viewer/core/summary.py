"""Dataset summary and quality control of a layer, without any UI (used by the
HTTP API's Summary panel)."""


def general_stats(sample, kind=None):
    """Observations with counts and UMIs of a sample. Cells layers also report the
    segmented cells, of which those without expression are kept on the map only."""
    adata = sample.adata
    mean_umis = 0.0
    if adata is not None and "total_counts" in adata.obs.columns:
        mean_umis = round(float(adata.obs["total_counts"].mean()), 1)
    area = sample.get_unit_area_um2()
    observations = int(adata.n_obs) if adata is not None else 0
    stats = dict(observations=observations,
                 mean_umis=mean_umis, umis_per_um2=round(float(mean_umis / area if area > 0 else 0), 1))
    if kind == "cells":
        segmented = int(len(sample.gdf))
        stats.update(segmented=segmented, not_sequenced=segmented - observations,
                     pct_sequenced=round(100 * observations / segmented, 1) if segmented else 0.0)
    return stats


def ensure_qc(adata):
    """Add scanpy's QC metrics (with mitochondrial genes) to `adata.obs` once."""
    import scanpy as sc
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    if "pct_counts_mt" not in adata.obs.columns:
        # scanpy's default top-gene proportions (50...500 genes) fail on smaller (targeted) panels
        top = [n for n in (50, 100, 200, 500) if n <= adata.n_vars] or None
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=top, inplace=True, log1p=False)
    return adata.obs
