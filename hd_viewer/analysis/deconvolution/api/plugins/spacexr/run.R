suppressPackageStartupMessages({library(Matrix); library(spacexr); library(jsonlite)})
if (!all(c('Reference', 'SpatialRNA', 'create.RCTD', 'run.RCTD') %in% getNamespaceExports('spacexr'))) {
    stop('This plugin needs the dmcable/spacexr GitHub API. Install it with remotes::install_github("dmcable/spacexr").')
}
folder <- commandArgs(trailingOnly=TRUE)[[1]]
spec <- fromJSON(file.path(folder, 'spec.json'))
p <- spec$params
set.seed(p$seed)

reference_counts <- as(readMM(file.path(folder, 'reference.mtx')), 'CsparseMatrix')
target_counts <- as(readMM(file.path(folder, 'target.mtx')), 'CsparseMatrix')
dimnames(reference_counts) <- list(spec$genes, spec$reference_ids)
dimnames(target_counts) <- list(spec$target_genes, spec$target_ids)
cell_types <- factor(spec$reference_labels, levels=spec$type_ids)
names(cell_types) <- spec$reference_ids
reference <- Reference(reference_counts, cell_types,
                       nUMI=colSums(reference_counts), min_UMI=p$reference_min_UMI,
                       n_max_cells=p$n_max_cells)
if (any(table(reference@cell_types) < p$CELL_MIN_INSTANCE)) {
    stop('Too few reference cells per type after UMI filtering. Check reference_min_UMI and CELL_MIN_INSTANCE.')
}

if (is.null(spec$coordinates)) {
    puck <- SpatialRNA(NULL, target_counts, use_fake_coords=TRUE)
} else {
    coords <- as.data.frame(spec$coordinates)
    names(coords) <- c('x', 'y')
    rownames(coords) <- spec$target_ids
    puck <- SpatialRNA(coords, target_counts)
}
# Report all filtered locations as unassigned rather than silently losing rows.
totals <- colSums(target_counts)
eligible <- totals >= p$UMI_min & totals <= p$UMI_max
if (!any(eligible)) {
    write_json(list(barcodes=character(), predictions=list(), genes=character()),
               file.path(folder, 'rctd.json'), auto_unbox=FALSE)
    quit(save='no', status=0)
}

cat('##progress 0.20 Selecting RCTD genes and estimating reference profiles\n')
fitted <- create.RCTD(
    puck, reference, max_cores=p$max_cores,
    gene_cutoff=p$gene_cutoff, fc_cutoff=p$fc_cutoff,
    gene_cutoff_reg=p$gene_cutoff_reg, fc_cutoff_reg=p$fc_cutoff_reg,
    UMI_min=p$UMI_min, UMI_max=p$UMI_max, counts_MIN=p$counts_MIN,
    UMI_min_sigma=p$UMI_min_sigma, CELL_MIN_INSTANCE=p$CELL_MIN_INSTANCE,
    MAX_MULTI_TYPES=p$MAX_MULTI_TYPES)

if (ncol(fitted@spatialRNA@counts) < 2L) {
    stop('RCTD needs at least two target observations after filtering for its dataset-level fit. Lower thresholds or use a larger layer.')
}
if (!any(fitted@spatialRNA@nUMI > p$UMI_min_sigma)) {
    stop('No target observations pass UMI_min_sigma. Lower the threshold or supply higher-depth counts.')
}

cat('##progress 0.40 Normalizing platform effects and fitting RCTD mixtures\n')
fitted <- run.RCTD(fitted, doublet_mode=p$doublet_mode)
results <- fitted@results
barcodes <- colnames(fitted@spatialRNA@counts)

# Preserve the native output, including doublet confidence / multi-mode details.
saveRDS(results, file.path(folder, '..', 'rctd_results.rds'))

predictions <- lapply(seq_along(barcodes), function(i) {
    barcode <- barcodes[[i]]
    if (p$doublet_mode == 'doublet') {
        row <- results$results_df[barcode, , drop=FALSE]
        list(spot_class=unbox(as.character(row$spot_class)),
             first_type=unbox(as.character(row$first_type)),
             second_type=unbox(as.character(row$second_type)),
             weights=as.numeric(results$weights_doublet[barcode, ]))
    } else if (p$doublet_mode == 'full') {
        list(types=colnames(results$weights), weights=as.numeric(results$weights[barcode, ]))
    } else {
        prediction <- results[[i]]
        list(types=as.character(prediction$cell_type_list), weights=as.numeric(prediction$sub_weights))
    }
})

genes <- unique(c(fitted@internal_vars$gene_list_reg, fitted@internal_vars$gene_list_bulk))
write_json(list(barcodes=barcodes, predictions=predictions, genes=genes),
           file.path(folder, 'rctd.json'), auto_unbox=FALSE, digits=NA, na='null')
           
cat('##progress 0.95 Exported RCTD predictions\n')
