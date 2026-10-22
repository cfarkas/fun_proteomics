#!/usr/bin/env Rscript
# Reformat the pinned DEP UbiLength data without changing measured LFQ values.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Usage: Rscript prepare_public_example.R UPSTREAM_DIR OUTPUT_DIR')
source <- normalizePath(args[[1]], mustWork=TRUE)
destination <- args[[2]]
if (dir.exists(destination) && length(list.files(destination, all.files=TRUE, no..=TRUE))) stop('Use a new or empty output directory')
dir.create(destination, recursive=TRUE, showWarnings=FALSE)
objects <- new.env(parent=emptyenv())
load(file.path(source, 'UbiLength.rda'), envir=objects)
load(file.path(source, 'UbiLength_ExpDesign.rda'), envir=objects)
data <- objects$UbiLength
design <- objects$UbiLength_ExpDesign
stopifnot(nrow(data)==3006L, nrow(design)==12L, !anyDuplicated(design$label))
columns <- paste0('LFQ.intensity.', design$label)
stopifnot(all(columns %in% names(data)))
flags <- c('Only.identified.by.site', 'Reverse', 'Potential.contaminant')
export <- data[, c('Protein.IDs', 'Gene.names', 'Protein.names', columns, flags)]
names(export)[1:3] <- c('Protein.Group', 'Genes', 'Protein.Names')
names(export)[4:15] <- design$label
samples <- data.frame(sample_id=design$label, condition=design$condition, replicate=design$replicate)
for (extension in c('csv', 'tsv')) {
  separator <- if (extension=='csv') ',' else '\t'
  write.table(export, file.path(destination, paste0('proteins.', extension)), sep=separator, row.names=FALSE, quote=TRUE, na='')
  write.table(samples, file.path(destination, paste0('samples.', extension)), sep=separator, row.names=FALSE, quote=TRUE, na='')
}
check <- read.csv(file.path(destination, 'proteins.csv'), check.names=FALSE)
stopifnot(identical(unname(as.matrix(check[, design$label])), unname(as.matrix(data[, columns]))))
cat('Exported 3,006 original rows and 12 samples. LFQ values and exclusion flags are unchanged.\n')
