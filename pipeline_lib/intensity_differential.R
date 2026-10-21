#!/usr/bin/env Rscript
# Missing intensities stay missing. Each named contrast has its own BH family.
suppressPackageStartupMessages({library(limma); library(jsonlite)})
job <- fromJSON(commandArgs(TRUE)[1])
samples <- read.csv(job$samples, check.names=FALSE, colClasses='character')
data <- read.csv(job$matrix, check.names=FALSE, colClasses='character')
if (!all(c('sample_id','condition','replicate') %in% names(samples))) stop('Samples require sample_id, condition, replicate')
if (anyDuplicated(samples$sample_id) || anyDuplicated(data[[1]])) stop('Duplicate sample or feature identifiers')
if (any(!nzchar(samples$condition)) || any(!nzchar(samples$replicate))) stop('Empty condition or replicate')
if (anyDuplicated(samples[c('condition','replicate')])) stop('Replicate IDs must be independent within each condition')
groups <- unique(samples$condition)
if (length(groups)<2 || any(table(samples$condition)<2)) stop('At least two independent replicates per condition and two conditions are required')
if (!all(samples$sample_id %in% names(data))) stop('Matrix is missing sample columns')
y <- as.matrix(data[,samples$sample_id,drop=FALSE]); storage.mode(y) <- 'numeric'
if (any(is.infinite(y)) || any(y<0,na.rm=TRUE)) stop('Intensities must be nonnegative finite values or missing')
y[y==0] <- NA_real_; y <- log2(y); rownames(y) <- data[[1]]
annotation <- data[,setdiff(names(data), samples$sample_id),drop=FALSE]; colnames(annotation)[1] <- 'FeatureID'
normalization <- if (is.null(job$normalization)) 'none' else job$normalization
offset <- rep(0,ncol(y))
if (normalization == 'median') {
 medians <- apply(y,2,median,na.rm=TRUE)
 if (any(!is.finite(medians))) stop('Median normalization requires observations in every sample')
 offset <- medians-median(medians); y <- sweep(y,2,offset)
} else if (normalization!='none') stop('Normalization must be none or median')
out <- job$outdir; dir.create(out, recursive=TRUE, showWarnings=FALSE)
write.csv(data.frame(sample_id=samples$sample_id,log2_offset=offset),file.path(out,'normalization.csv'),row.names=FALSE)
write.csv(cbind(annotation,as.data.frame(y,check.names=FALSE)),file.path(out,'log2_intensities.csv'),row.names=FALSE,na='')
design <- model.matrix(~0+factor(samples$condition,levels=groups)); colnames(design) <- groups
write.csv(data.frame(sample_id=samples$sample_id,design,check.names=FALSE),file.path(out,'design.csv'),row.names=FALSE)
mode <- if(is.null(job$comparisons)) 'all_pairwise' else job$comparisons
if (!mode %in% c('all_pairwise','control')) stop('Unknown comparison mode')
if (mode=='control' && (is.null(job$control)||!job$control %in% groups)) stop('Reference condition not found')
pairs <- if(mode=='control') lapply(setdiff(groups,job$control),function(g)c(job$control,g)) else combn(groups,2,simplify=FALSE)
labels <- make.unique(gsub('[^A-Za-z0-9_-]','_',groups),sep='_'); names(labels) <- groups
minimum <- if(is.null(job$min_replicates)) 2L else job$min_replicates
if(!is.finite(minimum)||minimum<2||minimum!=round(minimum)) stop('min_replicates must be an integer >= 2')
fdr <- if(is.null(job$fdr)) 0.05 else job$fdr
if(!is.finite(fdr)||fdr<=0||fdr>=1) stop('FDR must be between zero and one')
index <- list()
for(pair in pairs) {
 a <- samples$condition==pair[1]; b <- samples$condition==pair[2]
 na <- rowSums(is.finite(y[,a,drop=FALSE])); nb <- rowSums(is.finite(y[,b,drop=FALSE])); eligible <- na>=minimum & nb>=minimum
 tab <- annotation; tab$baseline <- rep(pair[1],nrow(y)); tab$comparison <- rep(pair[2],nrow(y))
 tab$logFC <- rowMeans(y[,b,drop=FALSE],na.rm=TRUE)-rowMeans(y[,a,drop=FALSE],na.rm=TRUE)
 tab$logFC[!is.finite(tab$logFC)] <- NA_real_; tab$PValue <- tab$FDR <- rep(NA_real_,nrow(y))
 tab$detected_reference <- na; tab$detected_comparison <- nb; tab$test_status <- rep('insufficient observations',nrow(y))
 if (sum(eligible)>=2) {
  # Pair-specific models tolerate missing values in unrelated conditions.
  paired <- a|b; d <- cbind(reference=as.integer(a[paired]),comparison=as.integer(b[paired]))
  fit <- lmFit(y[eligible,paired,drop=FALSE],d); fit <- contrasts.fit(fit,c(-1,1)); fit <- eBayes(fit,trend=TRUE)
  p <- as.numeric(fit$p.value[,1]); valid <- is.finite(p) & fit$df.residual>0; positions <- which(eligible)[valid]
  tab$logFC[positions] <- fit$coefficients[valid,1]; tab$PValue[positions] <- p[valid]
  tab$FDR[positions] <- p.adjust(p[valid],'BH'); tab$test_status[positions] <- 'tested'
 } else tab$test_status[eligible] <- 'insufficient features for variance moderation'
 tab$significant <- !is.na(tab$FDR)&tab$FDR<fdr; tab$method <- rep('limma-trend',nrow(tab)); tab$direction <- rep(paste(pair[2],'minus',pair[1]),nrow(tab))
 name <- paste0(labels[[pair[1]]],'_vs_',labels[[pair[2]]]); folder <- file.path(out,name); dir.create(folder,showWarnings=FALSE)
 write.csv(tab,file.path(folder,'all_results.csv'),row.names=FALSE,na=''); write.csv(tab[tab$significant,,drop=FALSE],file.path(folder,'significant_results.csv'),row.names=FALSE,na='')
 index[[length(index)+1]] <- data.frame(comparison=name,baseline=pair[1],comparison_condition=pair[2],method='limma-trend',reported_features=nrow(tab),tested_features=sum(!is.na(tab$PValue)),significant_features=sum(tab$significant),fdr=fdr,full_results=paste0(name,'/all_results.csv'),significant_results=paste0(name,'/significant_results.csv'))
}
write.csv(do.call(rbind,index),file.path(out,'results_index.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
