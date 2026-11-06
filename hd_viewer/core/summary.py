"""Dataset summary and quality control of a layer, without any UI (used by the
Bokeh Summary panel, board.layout.metadata, and the HTTP API)."""


def general_stats(sample):
    """Observation counts, sequencing rate and UMIs of a sample."""
    adata = sample.adata
    total = sample.n_geometries
    not_sequenced = sample.n_not_sequenced
    sequenced = total - not_sequenced
    mean_umis = 0.0
    if adata is not None and "total_counts" in adata.obs.columns:
        mean_umis = round(float(adata.obs["total_counts"].mean()), 1)
    area = sample.get_unit_area_um2()
    return dict(total=int(total), sequenced=int(sequenced), not_sequenced=int(not_sequenced),
                pct_sequenced=round(100 * sequenced / total, 1) if total > 0 else 0,
                mean_umis=mean_umis, umis_per_um2=round(float(mean_umis / area if area > 0 else 0), 1))


def ensure_qc(adata):
    """Add scanpy's QC metrics (with mitochondrial genes) to `adata.obs` once."""
    import scanpy as sc
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    if "pct_counts_mt" not in adata.obs.columns:
        # scanpy's default top-gene proportions (50...500 genes) fail on smaller (targeted) panels
        top = [n for n in (50, 100, 200, 500) if n <= adata.n_vars] or None
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=top, inplace=True, log1p=False)
    return adata.obs
