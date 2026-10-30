suppressPackageStartupMessages({library(Matrix); library(SingleCellExperiment); library(MuSiC); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
folder <- args[[1]]
spec <- fromJSON(file.path(folder, 'spec.json'))
p <- spec$params
reference <- readMM(file.path(folder, 'reference.mtx'))
target <- readMM(file.path(folder, 'target.mtx'))
rownames(reference) <- spec$genes
rownames(target) <- spec$target_genes
meta <- read.delim(file.path(folder, 'reference_obs.tsv'), row.names=1, check.names=FALSE, colClasses='character', na.strings=character())
colnames(reference) <- rownames(meta)
colnames(target) <- paste0('s', seq_len(ncol(target)) - 1L)
sce <- SingleCellExperiment(assays=list(counts=reference), colData=S4Vectors::DataFrame(meta))
parts <- list()
used.genes <- character()
for (start in seq(1L, ncol(target), by=p$batch_size)) {
    ix <- seq(start, min(ncol(target), start + p$batch_size - 1L))
    batch <- as.matrix(target[, ix, drop=FALSE])
    # MuSiC uses R subsetting without drop=FALSE internally. A duplicate keeps
    # singleton batches two-dimensional; each sample is fit independently.
    if (ncol(batch) == 1L) {
        batch <- cbind(batch, batch[, 1L])
        colnames(batch)[2L] <- '.dimension_guard'
    }
    fitted <- music_prop(bulk.mtx=batch, sc.sce=sce,
                        clusters=spec$label_column, samples=p$sample_column,
                        iter.max=p$iter_max, nu=p$nu, eps=p$eps, verbose=TRUE)
    used.genes <- union(used.genes, rownames(fitted$Weight.gene))
    value <- as.matrix(fitted$Est.prop.weighted)
    value <- value[colnames(target)[ix], , drop=FALSE]
    parts[[length(parts) + 1L]] <- value
    cat(sprintf('##progress %.4f MuSiC batch complete\n', max(ix)/ncol(target)))
}
fractions <- do.call(rbind, parts)
write.table(fractions, file.path(folder, 'fractions.tsv'), sep='\t', quote=TRUE, col.names=NA)
write_json(used.genes, file.path(folder, 'genes.json'), auto_unbox=FALSE)
