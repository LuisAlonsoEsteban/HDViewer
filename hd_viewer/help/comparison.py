"""Short explanations of the metrics returned by src.comparison.compare_samples.

"Spot" is a low-res geometry (sample 1); its "aggregate" is the area-weighted sum of the high-res
geometries (sample 2) that fall inside it.
"""

GLOBAL_HELP = {
    "tissue_area_sample1": "Total area of the sequenced low-res geometries (gaps between them do not count).",
    "tissue_area_sample2": "Total area of the sequenced high-res geometries (gaps between them do not count).",
    "tissue_intersection_area": "Area covered by sequenced geometries of both samples.",
    "tissue_dice": "2 x overlap / (area1 + area2). 1 = both cover the same area. Spots only cover ~1/4 of "
                   "the tissue, so spots vs. HD bins stay low (~0.35) even when perfectly aligned.",
    "n_sample1": "Sequenced low-res geometries.",
    "n_matched": "Low-res geometries that pass the coverage filter and have counts in both samples.",
    "n_sample2_matched": "High-res geometries used across the matched low-res geometries.",
    "n_genes": "Genes shared by both samples and detected in at least one of them.",
    "n_hvg": "Highly variable genes used by the HVG-based metrics.",
    "pseudobulk_pearson": "Correlation of each gene's overall expression in both samples (space ignored). "
                          "Low values point to gene naming, chemistry or probe-set differences.",
    "capture_ratio": "Total high-res counts / total low-res counts over the matched area: relative "
                     "sensitivity of the two technologies.",
    "total_counts_pearson": "Correlation of log total counts per spot: are RNA-rich spots rich in both "
                            "samples? Checks registration and tissue quality.",
    "pearson_hvg_median": "Median of the per-spot HVG correlation: main summary of whether matched spots "
                          "look the same.",
    "jensen_shannon_median": "Median per-spot Jensen-Shannon distance between expression compositions "
                             "(0 = identical, 1 = disjoint).",
    "gene_spatial_pearson_median": "Median, over genes detected in >=10% of spots in both samples, of the "
                                   "correlation of each gene across matched spots: same spatial maps?",
    "gene_spatial_pearson_median_hvg": "Same as gene_spatial_pearson_median, restricted to HVGs.",
    "morans_i_pearson": "Correlation of each HVG's spatial autocorrelation (Moran's I) between samples: "
                        "are the same genes spatially patterned?",
    "svg_top_jaccard": "Overlap (Jaccard) of the top spatially variable genes of each sample.",
    "local_metric_morans_i": "Moran's I of a per-spot metric map. High = disagreement is regional "
                             "(misregistration, damage, edges); ~0 = scattered noise.",
    "cluster_ari": "Adjusted Rand index between low-res clusters and the majority high-res cluster of "
                   "each spot. Works even if label names differ.",
    "cluster_match_rate": "Fraction of spots whose cluster equals the high-res majority cluster. Only "
                          "meaningful when both samples share label names.",
    "cluster_sample2_purity_mean": "Mean share of the dominant high-res cluster inside a spot. Low = spots "
                                   "mix populations, which limits any agreement.",
    "celltype_jensen_shannon_mean": "Mean divergence between deconvolved spot proportions and the "
                                    "high-res cell-type composition (0 = identical).",
    "celltype_spatial_pearson": "Per cell type, correlation of its proportion across spots in both "
                                "samples: does its spatial map agree?",
}

LOCAL_HELP = {
    "total_counts_sample1": "Total counts of the spot.",
    "total_counts_sample2": "Total counts of the aggregated high-res geometries inside the spot.",
    "log2_count_ratio": "log2(high-res counts / spot counts). Highlights regions where one technology "
                        "captures relatively more RNA.",
    "detected_genes_sample1": "Genes detected in the spot.",
    "detected_genes_sample2": "Genes detected in the high-res aggregate.",
    "pearson_hvg": "Correlation between the log-normalised HVG profiles of the spot and its high-res "
                   "aggregate. Main per-spot agreement score.",
    "jensen_shannon": "Jensen-Shannon distance between both expression compositions over all genes "
                      "(0 = identical, 1 = disjoint). Sensitive to depth and dropouts.",
    "coverage": "Fraction of the spot area covered by sequenced high-res geometries.",
    "n_sample2": "Number of high-res geometries inside the spot.",
    "distance_to_edge": "Distance from the spot to the tissue outline. Plot other metrics against it to "
                        "spot edge effects.",
    "cluster_sample2_majority": "Most abundant high-res cluster inside the spot.",
    "cluster_sample2_purity": "Share of the majority high-res cluster inside the spot.",
    "cluster_match": "Whether the spot's cluster equals the high-res majority cluster.",
    "celltype_jensen_shannon": "Divergence between deconvolved spot proportions and the high-res "
                               "cell-type composition (0 = identical).",
}

GENE_HELP = {
    "spatial_pearson": "Correlation of the gene across matched spots: does it have the same spatial map?",
    "capture_ratio": "High-res / low-res total counts of the gene. Outliers flag probe or chemistry biases.",
    "morans_i_sample1": "Spatial autocorrelation of the gene in the low-res sample (HVGs only).",
    "morans_i_sample2": "Spatial autocorrelation of the gene in the high-res aggregates (HVGs only).",
    "zero_fraction_sample1": "Fraction of matched spots where the gene is not detected (low-res).",
    "zero_fraction_sample2": "Fraction of matched spots where the gene is not detected (high-res).",
}


def global_help(metric):
    """Help of a flattened global metric name, e.g. 'local_metric_morans_i [pearson_hvg]'."""
    key, _, sub = metric.partition(" [")
    text = GLOBAL_HELP.get(key, "")
    if sub and key == "local_metric_morans_i":
        text = f"{text} Map: {sub.rstrip(']')}."
    elif sub and key == "celltype_spatial_pearson":
        text = f"{text} Cell type: {sub.rstrip(']')}."
    return text
