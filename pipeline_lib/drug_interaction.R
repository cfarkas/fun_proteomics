#!/usr/bin/env Rscript
# A replicated four-condition count contrast, not a viability synergy score.
suppressPackageStartupMessages({library(jsonlite);library(edgeR)})
script <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"drug_interaction_models.R"))

run_interactions <- function(job) {
 drugs <- job$optional_analyses$drugs
 opt <- drugs$interaction
 if (!isTRUE(opt$enabled)) return(invisible(NULL))
 if (!isTRUE(drugs$applied) || !identical(drugs$study_type, "synergy"))
  stop("Interaction analysis requires Drugs applied = Yes and a synergy study")
 methods <- if(is.null(opt$methods)) "edgeR" else unlist(opt$methods,use.names=FALSE)
 if(!length(methods)||anyDuplicated(methods)||any(!methods %in% c("edgeR","limma","monte_carlo")))
  stop("Choose one or more different interaction methods: edgeR, limma, monte_carlo")
 if("limma" %in% methods&&!requireNamespace("limma",quietly=TRUE))stop("limma is required for limma-voom interaction analysis")
 mc <- opt$monte_carlo
 integer_setting <- function(value,default,minimum,maximum,name) {
  if(is.null(value))value <- default
  if(!is.numeric(value)||length(value)!=1L||!is.finite(value)||value!=round(value)||value<minimum||value>maximum)
   stop(paste("Invalid Monte Carlo",name))
  as.integer(value)
 }
 if("monte_carlo" %in% methods) {
  mc <- list(simulations=integer_setting(mc$simulations,9999,99,1000000,"simulation count"),
             seed=integer_setting(mc$seed,1,0,2147483646,"seed"),
             workers=integer_setting(mc$workers,1,1,8,"workers"))
  available <- if(is.null(job$resources$cpus))1L else job$resources$cpus
  mc$effective_workers <- min(mc$workers,available,if(.Platform$OS.type=="windows")1L else 8L)
 }
 matching_confirmed <- isTRUE(opt$matched_design)
 interpretation <- if (matching_confirmed) "proteomic interaction" else "exploratory combination contrast"
 design_note <- if (matching_confirmed) {
  "Matching single-drug doses and treatment times were confirmed in the study settings."
 } else {
  "Matching single-drug doses and treatment times were not confirmed. Results are an exploratory combination contrast: dose or timing differences can contribute to the contrast, so it does not establish drug interaction in isolation. This annotation does not change the statistical test or its p-values."
 }
 roles <- c("control", "drug_a", "drug_b", "combination")
 conditions <- vapply(roles, function(role) opt[[role]], character(1))
 if (any(!nzchar(conditions)) || anyDuplicated(conditions)) stop("Select four distinct conditions for interaction analysis")
 cutoff <- if (is.null(opt$fdr)) 0.05 else opt$fdr
 if (!is.finite(cutoff) || cutoff <= 0 || cutoff > 1) stop("Invalid interaction FDR cutoff")
 parent <- file.path(job$outdir,"optional_analyses")
 dir.create(parent,recursive=TRUE,showWarnings=FALSE)
 target <- file.path(parent,"drug_interaction")
 base <- tempfile(".drug_interaction_",tmpdir=parent)
 dir.create(base)
 on.exit(unlink(base,recursive=TRUE),add=TRUE)
 summaries <- list()
 model_names <- c(edgeR="edgeR quasi-likelihood",limma="limma-voom empirical Bayes",
                  monte_carlo="negative-binomial null-model Monte Carlo")
 for (g in job$groups) {
  if (!identical(g$cohort, opt$cohort)) next
  slug <- if (g$mode == "DDA+") "ddaplus" else "dda"
  samples <- read.csv(file.path(g$directory, "samples.csv"), check.names=FALSE)
  all_samples <- samples$sample_id
  samples <- samples[samples$condition %in% conditions, , drop=FALSE]
  group <- factor(samples$condition, levels=conditions)
  if (any(table(group) < 2)) stop("Interaction analysis needs at least two biological replicates in each selected condition")
  if (any(vapply(split(samples$replicate, group), function(x) anyDuplicated(x)>0, logical(1))))
   stop("Interaction replicate numbers must be unique within each condition")
  design <- model.matrix(~0+group)
  colnames(design) <- roles
  contrast <- c(1, -1, -1, 1)
  null_design <- cbind(Intercept=1,drug_a=as.integer(group %in% conditions[c(2,4)]),drug_b=as.integer(group %in% conditions[c(3,4)]))
  kinds <- if (identical(job$analysis, "immunopeptidomics")) "peptide" else
   c("protein", if (!identical(job$optional_analyses$gene_counts,FALSE)) "gene", "peptide")
  for (kind in kinds) {
   folder <- file.path(g$cohort,slug,kind)
   common <- file.path(base,folder)
   dir.create(common,recursive=TRUE,showWarnings=FALSE)
   data <- read.csv(file.path(g$directory, paste0(kind, "_counts.csv")), check.names=FALSE)
   counts <- as.matrix(data[, samples$sample_id, drop=FALSE]); storage.mode(counts) <- "numeric"
   rownames(counts) <- data$FeatureID
   if (anyDuplicated(rownames(counts)) || any(!is.finite(counts)) || any(counts<0) || any(counts!=round(counts)))
    stop("Interaction counts must have unique identifiers and finite nonnegative integer values")
   keep <- rowSums(counts)>0
   write.csv(data.frame(sample_id=samples$sample_id,condition=samples$condition,design,check.names=FALSE),file.path(common,"design_matrix.csv"),row.names=FALSE)
   write.csv(data.frame(role=roles,condition=conditions,coefficient=contrast),file.path(common,"contrast.csv"),row.names=FALSE)
   write_interaction_matrix(counts,file.path(common,"raw_counts.csv"))
   normalized <- matrix(NA_real_,nrow(counts),ncol(counts),dimnames=dimnames(counts))
   normalization_status <- "insufficient evidence: no nonzero features or empty library"
   d <- NULL
   if(any(keep)&&all(colSums(counts)>0)) {
    d <- tryCatch(calcNormFactors(DGEList(counts[keep,,drop=FALSE],group=group),method="TMM"),
                  error=function(e){normalization_status<<-paste("insufficient evidence:",conditionMessage(e));NULL})
    if(!is.null(d)) {
     effective <- d$samples$lib.size*d$samples$norm.factors
     normalized <- cpm(counts,lib.size=effective)
     write.csv(data.frame(sample_id=colnames(d),library_size=d$samples$lib.size,TMM_factor=d$samples$norm.factors,
                         effective_library_size=effective,log_offset=log(effective)),file.path(common,"normalization_factors.csv"),row.names=FALSE)
    }
   }
   write_interaction_matrix(normalized,file.path(common,"model_normalized_CPM.csv"))
   for(method in methods) {
    cat("[DRUG INTERACTION]",g$cohort,g$mode,kind,method,"\n");flush.console()
    # Keep the existing edgeR location for older reports and links.
    relative <- if(method=="edgeR")folder else file.path(folder,method)
    out <- file.path(base,relative)
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    status <- normalization_status
    fitted <- NULL
    if(!is.null(d)) {
     status <- "tested"
     fitted <- tryCatch(switch(method,
      edgeR=interaction_edger(d,design,contrast,out),
      limma=interaction_limma(d,design,contrast,out),
      monte_carlo=interaction_monte_carlo(d,design,null_design,contrast,out,mc$simulations,
       interaction_seed(mc$seed,paste(g$cohort,g$mode,kind,sep="/")),mc$effective_workers)),
      error=function(e){status<<-paste("insufficient evidence:",conditionMessage(e));NULL})
    }
    if(is.null(fitted)) {
     writeLines(status,file.path(out,"INSUFFICIENT_EVIDENCE.txt"))
     tab <- data.frame(interaction_log2FC=rep(NA_real_,nrow(counts)),PValue=rep(NA_real_,nrow(counts)),row.names=rownames(counts))
     tab$test_status <- rep(status,nrow(tab))
    } else {
     # Reinsert all-zero rows without inventing test results.
     tab <- fitted[match(rownames(counts),rownames(fitted)),,drop=FALSE]
     rownames(tab) <- rownames(counts)
     if(is.null(tab$test_status))tab$test_status <- ifelse(is.finite(tab$PValue),"tested","insufficient evidence: nonfinite model statistic")
     tab$test_status[!keep] <- "insufficient evidence: all_zero"
    }
    valid <- is.finite(tab$PValue)&tab$PValue>=0&tab$PValue<=1
    tab$PValue[!valid] <- NA_real_
    tab$test_status[is.na(tab$test_status)] <- "insufficient evidence: nonfinite model statistic"
    tab$FDR <- p.adjust(tab$PValue,method="BH")
    tab$significant <- !is.na(tab$FDR)&tab$FDR<cutoff
    tab$direction <- ifelse(!is.finite(tab$interaction_log2FC),"not_estimable",
     ifelse(tab$interaction_log2FC>0,"above_log_additive_expectation",
      ifelse(tab$interaction_log2FC<0,"below_log_additive_expectation","at_log_additive_expectation")))
    tab$method <- rep(method,nrow(tab));tab$model <- rep(model_names[[method]],nrow(tab))
    tab$dose_time_matching_confirmed <- rep(matching_confirmed,nrow(tab))
    tab$analysis_interpretation <- rep(interpretation,nrow(tab))
    for(i in seq_along(roles)) {
     selected <- group==conditions[i]
     tab[[paste0("mean_count_",roles[i])]] <- rowMeans(counts[,selected,drop=FALSE])
     tab[[paste0("mean_CPM_",roles[i])]] <- rowMeans(normalized[,selected,drop=FALSE])
     tab[[paste0("detected_samples_",roles[i])]] <- rowSums(counts[,selected,drop=FALSE]>0)
     tab[[paste0("total_samples_",roles[i])]] <- rep(sum(selected),nrow(tab))
    }
    write.csv(data.frame(FeatureID=as.character(rownames(counts)),tested=valid,reason=tab$test_status),file.path(out,"feature_filter.csv"),row.names=FALSE)
    annotation <- data[match(rownames(tab),data$FeatureID),setdiff(names(data),all_samples),drop=FALSE]
    tab <- cbind(annotation,tab)
    tab <- tab[order(tab$PValue,na.last=TRUE),,drop=FALSE]
    write.csv(tab,file.path(out,"interaction.csv"),row.names=FALSE,na="")
    write.csv(tab[tab$significant,,drop=FALSE],file.path(out,"significant.csv"),row.names=FALSE,na="")
    if(!any(valid)&&status=="tested")status <- "insufficient evidence: no estimable tests"
    if(any(valid)&&any(keep&!valid))status <- "tested with some unestimable features"
    summaries[[length(summaries)+1L]] <- data.frame(cohort=g$cohort,mode=g$mode,feature_type=kind,method=method,
     input_features=nrow(counts),tested_features=sum(valid),significant_features=sum(tab$significant),status=status,
     results_file=file.path(relative,"interaction.csv"),dose_time_matching_confirmed=matching_confirmed,
     analysis_interpretation=interpretation,MC_minimum_pvalue=if(method=="monte_carlo")1/(mc$simulations+1) else NA_real_)
   }
  }
 }
 if(!length(summaries))stop("No analysis cohort matches the interaction settings")
 write.csv(do.call(rbind,summaries),file.path(base,"analysis_summary.csv"),row.names=FALSE,na="")
 opt$methods <- as.list(methods)
 if("monte_carlo" %in% methods)opt$monte_carlo <- mc else opt$monte_carlo <- NULL
 settings <- c(opt,list(dose_time_matching_confirmed=matching_confirmed,analysis_interpretation=interpretation))
 write_json(settings,file.path(base,"settings.json"),auto_unbox=TRUE,pretty=TRUE)
 writeLines(capture.output(sessionInfo()),file.path(base,"R_session.txt"))
 writeLines(c("# Drug-combination interaction", "",
  "All selected methods test combination - drug A - drug B + control on a log2 scale, using the same TMM-normalized spectral-count libraries and four-condition design. The null is multiplicative on the original abundance scale. Independent biological replicates are assumed; batch and paired-subject effects are not inferred.","",
  "edgeR uses a negative-binomial quasi-likelihood model. limma-voom uses precision-weighted log2 CPM and empirical Bayes variance moderation, with moderated t-tests and 95% confidence intervals. limma's small-count transformation can yield different effect estimates from the count GLMs; no counts are imputed. Empirical Bayes moderation is not a fully Bayesian interaction posterior.","",
  "Monte Carlo fits a negative-binomial model with control and both single-drug main effects but no interaction. It simulates counts from those fitted null means, retaining fitted dispersions and effective library offsets. Both null and full coefficient models are refitted to every simulated dataset. The test statistic is the likelihood-ratio deviance; its two-sided departure includes either direction of interaction. PValue = (1 + number of simulated statistics at least as large as observed) / (1 + simulations). This is a plug-in parametric bootstrap, conditional on estimated dispersions and normalization, not a permutation test or a Bayesian posterior. Its calibration is approximate, especially with few counts or replicates.","",
  "Monte Carlo tables include exceedance counts, simulation totals, the derived family seed, an approximate simulation standard error and a 95% exact binomial interval for the conditional tail probability. Those intervals measure Monte Carlo uncertainty, not uncertainty in biological effect size or fitted nuisance parameters. Minimum attainable p-value is 1/(simulations+1); 9,999 simulations gives 0.0001 and may limit BH discoveries in large feature families. More simulations improve p-value resolution, not biological sample size. Independent reproducible RNG streams give the same result with different worker counts on the same input and software. The base seed and effective worker cap are saved in settings.json; null means, dispersions, design and offsets are exported.","",design_note,"",
  "Positive interaction_log2FC means above the log-additive expectation; negative means below. Direction alone does not mean beneficial synergy or antagonism. The contrast describes proteomic responses at the measured doses and time, not a Bliss, Loewe, ZIP or cell-viability synergy score. HLA analyses test peptide counts only and do not measure source-protein abundance.","",
  "Each method, cohort, search mode and feature type has a separate test family with Benjamini-Hochberg FDR. Do not select the smallest p-value across methods as a combined test. All nonzero features enter each model. Low coverage does not fail the study: full tables retain nonsignificant and all-zero features; untestable rows have an explicit test_status, blank p-values/FDR and significant=false. Untestable rows are not evidence of no effect. A model that cannot be fitted does not prevent other selected methods from running. No pseudocounts or imputation are added to the count-model input. The analysis runs independently of the DEPs switch.","",
  "analysis_summary.csv lists each method's result path. Existing edgeR results remain at <cohort>/<mode>/<feature_type>/interaction.csv; limma and monte_carlo use corresponding subdirectories. Each has significant.csv and feature_filter.csv. Shared raw counts, design, contrast and TMM normalization are in the feature-type directory. limma saves voom log2 CPM and weights. Only selected methods appear in a newly generated report; a completed report replaces previous interaction outputs."),file.path(base,"README.md"))
 # Publish only once all selected methods have produced their reports.
 if(dir.exists(target))unlink(target,recursive=TRUE)
 if(!file.rename(base,target))stop("Could not publish drug-interaction results")
 invisible(target)
}
run_interactions(fromJSON(commandArgs(TRUE)[1],simplifyVector=FALSE))
