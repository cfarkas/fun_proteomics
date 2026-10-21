#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(jsonlite);library(edgeR);library(limma)})
if(packageVersion("edgeR") < "4.4.2")stop("edgeR >= 4.4.2 is required for sparse-count quasi-likelihood modeling")
if(packageVersion("limma") < "3.62.2")stop("limma >= 3.62.2 is required")
args<-commandArgs(TRUE);job<-fromJSON(args[1],simplifyVector=FALSE);opt<-job$differential
methods<-unlist(opt$methods);if(is.null(methods))methods<-c("edgeR","limma")
if(isTRUE(opt$enabled)&&"edgeRun" %in% methods&&!requireNamespace("edgeRun",quietly=TRUE))stop("edgeRun is not installed; install dependencies before running")
upper<-if(is.null(opt$edgeRun_upper))100000 else as.integer(opt$edgeRun_upper)
if(!is.finite(upper)||upper<100000||upper%%100!=0)stop("edgeRun upper must be at least 100000 and a multiple of 100")
workers<-if(is.null(opt$workers))1L else as.integer(opt$workers)
if(requireNamespace("data.table",quietly=TRUE))data.table::setDTthreads(1)
fdr<-if(is.null(opt$fdr))0.05 else opt$fdr
control<-if(is.null(opt$control))"CTRL" else opt$control
allpairs<-identical(opt$comparisons,"all_pairwise")||identical(control,"")
write_matrix<-function(m,path){write.csv(data.frame(FeatureID=as.character(rownames(m)),m,check.names=FALSE),path,row.names=FALSE,na="")}
write_result<-function(tab,annotation,path,method,contrast){
 valid<-!is.na(tab$PValue)
 if(any(!is.finite(tab$PValue[valid]))||any(tab$PValue[valid] < -1e-8 | tab$PValue[valid] > 1+1e-8))stop("Invalid test p-values")
 if(!"test_status" %in% names(tab))tab$test_status<-ifelse(valid,"tested","insufficient evidence")
 tab$PValue<-pmin(1,pmax(0,tab$PValue))
 tab$FeatureID<-rownames(tab);tab$adjusted_p_value<-p.adjust(tab$PValue,method="BH");tab$significant<-!is.na(tab$adjusted_p_value)&tab$adjusted_p_value<fdr
 tab$method<-rep(method,nrow(tab));tab$contrast<-rep(contrast,nrow(tab))
 if(isTRUE(job$retain_untested_rows)){
  missing<-setdiff(annotation$FeatureID,tab$FeatureID)
  if(length(missing)){
   extra<-tab[rep(NA_integer_,length(missing)),,drop=FALSE];extra$FeatureID<-missing
   extra$test_status<-ifelse(rowSums(counts[missing,,drop=FALSE])==0,"not tested: all-zero counts","not tested: low-count filter")
   for(column in intersect(c("baseline","comparison","direction","count_unit","normalization","model","FeatureType"),names(tab)))extra[[column]]<-tab[[column]][1]
   extra$significant<-FALSE;extra$method<-method;extra$contrast<-contrast;tab<-rbind(tab,extra)
  }
 }
 ann<-annotation[match(tab$FeatureID,annotation$FeatureID),,drop=FALSE]
 tab<-cbind(ann,tab[,setdiff(names(tab),"FeatureID"),drop=FALSE]);tab<-tab[order(tab$PValue,na.last=TRUE),]
 write.csv(tab,path,row.names=FALSE,na="");return(tab)
}
summary<-list();comparison_rows<-list();analysis_rows<-list();method_notes<-list()
filter_counts<-isTRUE(opt$enabled)&&isTRUE(opt$filter_low_counts)
min_count<-if(is.null(opt$min_count))10 else opt$min_count
min_total<-if(is.null(opt$min_total_count))15 else opt$min_total_count
if(filter_counts&&(!is.finite(min_count)||min_count<1||!is.finite(min_total)||min_total<1))stop("Count filtering thresholds must be positive when low-count filtering is enabled")
mod_min_count<-if(is.null(opt$modification_min_count))1 else opt$modification_min_count
mod_min_total<-if(is.null(opt$modification_min_total_count))3 else opt$modification_min_total_count
if(filter_counts&&isTRUE(opt$modifications)&&(!is.finite(mod_min_count)||mod_min_count<1||!is.finite(mod_min_total)||mod_min_total<1))stop("Modification filtering thresholds must be positive when modification filtering is enabled")
categories<-c(protein="proteins",gene="genes",peptide="peptides",modification_gene="modifications_by_gene",modification_protein="modifications_by_protein")
custom_categories<-unlist(job$feature_categories)
if(length(custom_categories))categories<-custom_categories
sparse_kind<-function(kind)startsWith(kind,"modification")||kind %in% unlist(job$sparse_kinds)
for(g in job$groups){
 cat("[FINAL RESULTS]",g$cohort,g$mode,"normalization\n")
 folder<-g$directory;samples<-read.csv(file.path(folder,"samples.csv"),check.names=FALSE)
 group<-factor(samples$condition,levels=unique(samples$condition))
 slug<-if(g$mode=="DDA+")"ddaplus" else "dda"
 normalized<-file.path(job$outdir,"matrices",g$cohort,slug)
 zdir<-file.path(dirname(job$outdir),"zscore",g$cohort,slug)
 if(!is.null(job$zscore_directory))zdir<-file.path(job$zscore_directory,g$cohort,slug)
 dir.create(normalized,recursive=TRUE,showWarnings=FALSE);dir.create(zdir,recursive=TRUE,showWarnings=FALSE)
 matrices<-list();annotations<-list()
 for(kind in if(length(custom_categories))names(custom_categories) else c("protein","peptide","modification","gene","modification_gene","modification_protein")){
  data<-read.csv(file.path(folder,paste0(kind,"_counts.csv")),check.names=FALSE)
  counts<-as.matrix(data[,samples$sample_id,drop=FALSE]);storage.mode(counts)<-"numeric";rownames(counts)<-data$FeatureID
  if(any(!is.finite(counts))||any(counts<0)||any(counts!=round(counts)))stop("Counts must be finite nonnegative integers")
  matrices[[kind]]<-counts;annotations[[kind]]<-data[,setdiff(names(data),samples$sample_id),drop=FALSE]
  if(nrow(counts)==0){write_matrix(counts,file.path(normalized,paste0(kind,"_normalized_cpm.csv")));write_matrix(counts,file.path(normalized,paste0(kind,"_log2_cpm.csv")));write_matrix(counts,file.path(zdir,paste0(kind,"_zscore.csv")));next}
  positive<-colSums(counts)>0
  cpmat<-logmat<-matrix(NA_real_,nrow(counts),ncol(counts),dimnames=dimnames(counts))
  if(any(positive)){
   d<-DGEList(counts[,positive,drop=FALSE]);d<-calcNormFactors(d,method=if(sparse_kind(kind))"TMMwsp" else "TMM")
   cpmat[,positive]<-cpm(d);logmat[,positive]<-cpm(d,log=TRUE,prior.count=2)
   write.csv(data.frame(sample_id=colnames(d),library_size=d$samples$lib.size,TMM_factor=d$samples$norm.factors),file.path(folder,paste0(kind,"_normalization_factors.csv")),row.names=FALSE)
  }
  z<-logmat
  for(i in seq_len(nrow(z))){good<-is.finite(logmat[i,]);v<-logmat[i,good];z[i,good]<-if(length(v)<2||diff(range(v))<=1e-12)0 else (v-mean(v))/sd(v)}
  write_matrix(cpmat,file.path(normalized,paste0(kind,"_normalized_cpm.csv")));write_matrix(logmat,file.path(normalized,paste0(kind,"_log2_cpm.csv")));write_matrix(z,file.path(zdir,paste0(kind,"_zscore.csv")))

 }

 if(!isTRUE(opt$enabled))next
 if(nlevels(group)<2||any(table(group)<2))stop("Differential analysis requires at least two conditions and two independent replicates per condition")
 design<-model.matrix(~0+group);colnames(design)<-levels(group)
 levels<-levels(group)
 pairs<-if(!allpairs&&control %in% levels)lapply(setdiff(levels,control),function(x)c(control,x)) else combn(levels,2,simplify=FALSE)
 labels<-vapply(levels,function(x){value<-opt$condition_labels[[x]];if(is.null(value)||!nzchar(value))x else value},character(1))
 labels<-make.unique(gsub("[^A-Za-z0-9_-]","_",labels),sep="_");names(labels)<-levels
 contrasts<-sapply(pairs,function(pair){v<-rep(0,length(levels));v[match(pair[1],levels)]<--1;v[match(pair[2],levels)]<-1;v})
 if(is.null(dim(contrasts)))contrasts<-matrix(contrasts,ncol=1)
 rownames(contrasts)<-levels;colnames(contrasts)<-make.unique(vapply(pairs,function(pair)paste0(labels[[pair[1]]],"_vs_",labels[[pair[2]]]),character(1)),sep="_")
 for(i in seq_along(pairs))comparison_rows[[length(comparison_rows)+1]]<-data.frame(cohort=g$cohort,mode=g$mode,contrast=colnames(contrasts)[i],baseline=pairs[[i]][1],comparison=pairs[[i]][2],direction=paste(pairs[[i]][2],"minus",pairs[[i]][1]))
 # Protein comparisons are primary. Gene comparisons are explicitly optional.
 gene_counts<-!identical(job$optional_analyses$gene_counts,FALSE)
 test_kinds<-if(identical(job$analysis,"immunopeptidomics"))"peptide" else c("protein","peptide",if(isTRUE(opt$modifications))"modification_protein",if(gene_counts)c("gene",if(isTRUE(opt$modifications))"modification_gene"))
 if(length(custom_categories))test_kinds<-names(custom_categories)
 for(kind in test_kinds){
 counts<-matrices[[kind]];annotation<-annotations[[kind]]
 category<-categories[[kind]];sparse<-sparse_kind(kind)
 norm_method<-if(sparse)"TMMwsp" else "TMM"
 diag<-file.path(folder,"models",category);dir.create(diag,recursive=TRUE,showWarnings=FALSE)
 out<-file.path(job$outdir,"differential")
 result_path<-function(method,name){folder<-file.path(out,method,g$cohort,slug);dir.create(folder,recursive=TRUE,showWarnings=FALSE);file.path(folder,paste0(method,"_",name,"_",category,".csv"))}
 record_status<-function(tested,status){analysis_rows[[length(analysis_rows)+1]]<<-data.frame(Cohort=g$cohort,Mode=g$mode,FeatureType=kind,input_features=nrow(counts),tested_features=tested,zero_fraction=if(length(counts))mean(counts==0) else NA_real_,normalization=norm_method,filter_low_counts=filter_counts,min_count=if(!filter_counts)NA_real_ else if(sparse)mod_min_count else min_count,min_total_count=if(!filter_counts)NA_real_ else if(sparse)mod_min_total else min_total,status=status)}
 unavailable_result<-function(reason,selected_methods=unique(c("edgeR",methods)),only_contrast=NULL){
  status<-paste("insufficient evidence:",reason)
  writeLines(status,file.path(diag,"INSUFFICIENT_EVIDENCE.txt"))
  if(!file.exists(file.path(diag,"counts_model_normalized_cpm.csv")))file.copy(file.path(normalized,paste0(kind,"_normalized_cpm.csv")),file.path(diag,"counts_model_normalized_cpm.csv"),overwrite=TRUE)
  tab<-data.frame(logFC=rep(NA_real_,nrow(counts)),PValue=rep(NA_real_,nrow(counts)),test_status=rep(status,nrow(counts)),row.names=rownames(counts))
  for(method in selected_methods){
   if(is.null(only_contrast)&&method %in% c("edgeR","limma"))write_result(tab,annotation,result_path(method,if(method=="edgeR")"GLM_omnibus" else "omnibus"),method,"any_condition_difference")
   if(method %in% methods)for(i in seq_along(pairs))if(is.null(only_contrast)||colnames(contrasts)[i]==only_contrast)write_result(tab,annotation,result_path(method,colnames(contrasts)[i]),method,colnames(contrasts)[i])
  }
 }
 if(nrow(counts)==0||any(colSums(counts)==0)){
  unavailable_result("no features or an empty count library")
  write.csv(data.frame(FeatureID=as.character(rownames(counts)),tested=rep(FALSE,nrow(counts)),reason=rep("insufficient_evidence",nrow(counts))),file.path(diag,"counts_test_filter.csv"),row.names=FALSE)
  record_status(0,"insufficient evidence: no features or empty library");next
 }
 d<-DGEList(counts,group=group);if(sparse)d<-calcNormFactors(d,method=norm_method);nonzero<-rowSums(counts)>0;keep<-nonzero;if(filter_counts)keep<-nonzero & filterByExpr(d,design=design,min.count=if(sparse)mod_min_count else min_count,min.total.count=if(sparse)mod_min_total else min_total);write.csv(data.frame(FeatureID=rownames(d),tested=keep,reason=ifelse(!nonzero,"all_zero",ifelse(keep,"retained","low_count"))),file.path(diag,"counts_test_filter.csv"),row.names=FALSE)
 if(sum(keep)==0){unavailable_result("no nonzero features after filtering");record_status(0,"insufficient evidence after filtering");next}
 d<-d[keep,,keep.lib.sizes=sparse];if(!sparse)d<-calcNormFactors(d,method=norm_method)
 model_error<-NULL
 fit<-tryCatch({d<-estimateDisp(d,design,robust=TRUE);if(sparse)glmQLFit(d,design,robust=TRUE,legacy=FALSE) else glmQLFit(d,design,robust=TRUE,legacy=TRUE)},error=function(e){model_error<<-conditionMessage(e);NULL})
 if(is.null(fit)){
  unavailable_result(model_error,intersect(unique(c("edgeR",methods)),c("edgeR","edgeRun")))
  record_status(0,paste("insufficient evidence for edgeR:",model_error))
 }else record_status(sum(keep),"tested")
 write_matrix(cpm(d),file.path(diag,"counts_model_normalized_cpm.csv"))
 write_matrix(cpm(d,log=TRUE,prior.count=2),file.path(diag,"counts_model_log2_cpm.csv"))
 write.csv(data.frame(sample_id=colnames(d),library_size=d$samples$lib.size,TMM_factor=d$samples$norm.factors),file.path(diag,"counts_model_normalization_factors.csv"),row.names=FALSE)
 write.csv(data.frame(sample_id=samples$sample_id,design,check.names=FALSE),file.path(diag,"design_matrix.csv"),row.names=FALSE)
 write.csv(data.frame(condition=levels,contrasts,check.names=FALSE),file.path(diag,"contrast_matrix.csv"),row.names=FALSE)
 # GLM in every experiment, including two-condition studies. Overall test spans K-1 independent contrasts.
 omni<-diag(length(levels))[,seq.int(2,length(levels)),drop=FALSE];omni[1,]<--1
 if(!is.null(fit)){
  overall<-topTags(glmQLFTest(fit,contrast=omni),n=Inf,sort.by="none")$table
  write_result(overall,annotation,result_path("edgeR","GLM_omnibus"),"edgeR_GLM","any_condition_difference")
 }
 vfit<-NULL
 if("limma" %in% methods){
  v<-tryCatch(voom(d,design,plot=FALSE),error=function(e){
   reason<-paste("limma-voom not estimable:",conditionMessage(e));warning(reason);unavailable_result(reason,"limma")
   method_notes[[length(method_notes)+1]]<<-data.frame(Cohort=g$cohort,Mode=g$mode,FeatureType=kind,method="limma",status=reason)
   path<-result_path("limma","NOT_ESTIMABLE");writeLines(reason,sub("\\.csv$",".txt",path));NULL
  })
  if(!is.null(v)){
  vfit<-tryCatch({lf<-lmFit(v,design);eBayes(contrasts.fit(lf,contrasts),robust=TRUE)},error=function(e){unavailable_result(conditionMessage(e),"limma");method_notes[[length(method_notes)+1]]<<-data.frame(Cohort=g$cohort,Mode=g$mode,FeatureType=kind,method="limma",status=paste("insufficient evidence:",conditionMessage(e)));NULL})
  if(!is.null(vfit)){
  if(is.null(fit)){analysis_rows[[length(analysis_rows)]]$tested_features<-sum(keep);analysis_rows[[length(analysis_rows)]]$status<-"limma tested; edgeR insufficient evidence"}
  dimnames(v$weights)<-dimnames(v$E);write_matrix(v$E,file.path(diag,"limma_voom_log2_cpm.csv"));write_matrix(v$weights,file.path(diag,"limma_voom_weights.csv"))
  omfit<-eBayes(contrasts.fit(lf,omni),robust=TRUE);ot<-topTable(omfit,coef=seq_len(ncol(omni)),number=Inf,sort.by="none");ot$PValue<-ot$P.Value
  write_result(ot,annotation,result_path("limma","omnibus"),"limma","any_condition_difference")
  }
  }
 }
 for(i in seq_along(pairs)){
  name<-colnames(contrasts)[i]
  for(method in methods){
   if(method=="limma"&&is.null(vfit))next
   if(method %in% c("edgeR","edgeRun")&&is.null(fit))next
   cat("[DIFFERENTIAL]",g$cohort,g$mode,category,method,name,"\n");flush.console()
   if(method=="edgeR")tab<-topTags(glmQLFTest(fit,contrast=contrasts[,i]),n=Inf,sort.by="none")$table
   else if(method=="limma"){tab<-topTable(vfit,coef=i,number=Inf,sort.by="none");tab$PValue<-tab$P.Value}
   else if(method=="edgeRun"){
    # An unconditional two-group test, never labeled as a GLM.
    chunks<-split(seq_len(nrow(d)),ceiling(seq_len(nrow(d))/64))
    pieces<-parallel::mclapply(chunks,function(ix){dd<-d[ix,,keep.lib.sizes=TRUE];edgeRun::UCexactTest(dd,pair=pairs[[i]],dispersion=d$tagwise.dispersion[ix],upper=upper,by=100)$table},mc.cores=workers)
    if(any(vapply(pieces,inherits,logical(1),"try-error"))){
     reason<-paste("edgeRun could not estimate this contrast:",paste(unlist(pieces[vapply(pieces,inherits,logical(1),"try-error")]),collapse="; "))
     unavailable_result(reason,"edgeRun",name);method_notes[[length(method_notes)+1]]<-data.frame(Cohort=g$cohort,Mode=g$mode,FeatureType=kind,method="edgeRun",status=reason);next
    }
    tab<-do.call(rbind,unname(pieces));if(!all(rownames(d) %in% rownames(tab)))stop("edgeRun feature identifiers were lost");tab<-tab[rownames(d),,drop=FALSE]
   }else stop("Unknown differential method")
   tab$FeatureType<-kind;tab$count_unit<-if(!is.null(job$count_units[[kind]]))job$count_units[[kind]] else if(startsWith(kind,"modification_"))"distinct modified peptide sequences" else "spectral counts";tab$normalization<-norm_method;tab$model<-if(method=="edgeR")if(sparse)"edgeR v4 small-count quasi-likelihood GLM" else "edgeR quasi-likelihood GLM" else if(method=="edgeRun")"unconditional exact test" else "limma-voom"
   tab$baseline<-pairs[[i]][1];tab$comparison<-pairs[[i]][2];tab$direction<-paste(pairs[[i]][2],"minus",pairs[[i]][1])
   result<-write_result(tab,annotation,result_path(method,name),method,name)
   result$Cohort<-g$cohort;result$Mode<-g$mode;result$FeatureType<-kind
   result$baseline<-pairs[[i]][1];result$comparison<-pairs[[i]][2];result$direction<-paste(pairs[[i]][2],"minus",pairs[[i]][1])
   for(column in c("Protein","Gene","Peptide","Assigned Modifications","Modification"))if(!column %in% names(result))result[[column]]<-""
   summary[[length(summary)+1]]<-result[,c("Cohort","Mode","FeatureType","FeatureID","Protein","Gene","Peptide","Assigned Modifications","Modification","count_unit","normalization","model","method","contrast","baseline","comparison","direction","logFC","PValue","adjusted_p_value","significant")]
  }
 }
 }
}
if(length(summary)){
 results<-do.call(rbind,summary)
 for(method in unique(results$method)){
  selected<-results[results$method==method & results$significant,]
  write.csv(selected,file.path(job$outdir,"differential",method,"significant_results.csv"),row.names=FALSE)
 }

}
if(length(comparison_rows))write.csv(do.call(rbind,comparison_rows),file.path(job$outdir,"comparisons.csv"),row.names=FALSE)
if(length(method_notes))write.csv(do.call(rbind,method_notes),file.path(job$outdir,"method_notes.csv"),row.names=FALSE)
if(length(analysis_rows))write.csv(do.call(rbind,analysis_rows),file.path(job$outdir,"analysis_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(job$outdir,"R_session.txt"))
