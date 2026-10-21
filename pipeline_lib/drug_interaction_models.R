# Model helpers for the same combination - A - B + control contrast.
write_interaction_matrix <- function(x, path) {
 write.csv(data.frame(FeatureID=as.character(rownames(x)),x,check.names=FALSE),path,row.names=FALSE,na="")
}

interaction_edger <- function(d, design, contrast, out) {
 d <- estimateDisp(d,design,robust=TRUE)
 fit <- glmQLFit(d,design,robust=TRUE,legacy=TRUE)
 tab <- topTags(glmQLFTest(fit,contrast=contrast),n=Inf,sort.by="none")$table
 names(tab)[names(tab)=="logFC"] <- "interaction_log2FC"
 tab
}

interaction_limma <- function(d, design, contrast, out) {
 v <- limma::voom(d,design,plot=FALSE)
 cm <- matrix(contrast,ncol=1,dimnames=list(colnames(design),"interaction"))
 fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v,design),cm),robust=TRUE)
 tab <- limma::topTable(fit,coef=1,number=Inf,sort.by="none",confint=TRUE)
 names(tab)[names(tab)=="logFC"] <- "interaction_log2FC"
 names(tab)[names(tab)=="P.Value"] <- "PValue"
 names(tab)[names(tab)=="CI.L"] <- "interaction_CI95_lower"
 names(tab)[names(tab)=="CI.R"] <- "interaction_CI95_upper"
 tab$adj.P.Val <- NULL
 dimnames(v$weights) <- dimnames(v$E)
 write_interaction_matrix(v$E,file.path(out,"voom_log2_CPM.csv"))
 write_interaction_matrix(v$weights,file.path(out,"voom_weights.csv"))
 tab
}

interaction_seed <- function(seed, label) {
 # Stable per-family seeds keep subsets and worker scheduling reproducible.
 value <- as.double(seed)
 for(byte in utf8ToInt(label))value <- (value*131+byte) %% 2147483646
 as.integer(value)
}

interaction_monte_carlo <- function(d, design, null_design, contrast, out, simulations, seed, workers=1L) {
 d <- estimateDisp(d,design,robust=TRUE)
 dispersion <- d$tagwise.dispersion
 if(any(!is.finite(dispersion)|dispersion<0))stop("Invalid dispersion estimates for Monte Carlo simulation")
 fit <- glmFit(d,design,dispersion=dispersion,prior.count=0)
 null <- glmFit(d,null_design,dispersion=dispersion,prior.count=0)
 observed <- pmax(0,null$deviance-fit$deviance)
 if(any(!is.finite(observed))||any(!is.finite(null$fitted.values)))stop("Nonfinite fitted null model")
 offset <- getOffset(d)
 write.csv(data.frame(sample_id=colnames(d),null_design,check.names=FALSE),file.path(out,"null_design_matrix.csv"),row.names=FALSE)
 write_interaction_matrix(null$fitted.values,file.path(out,"null_fitted_counts.csv"))
 write.csv(data.frame(FeatureID=rownames(d),dispersion=dispersion,observed_LR=observed),file.path(out,"null_model_parameters.csv"),row.names=FALSE)

 old_kind <- RNGkind()
 had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
 if(had_seed)old_seed <- get(".Random.seed",envir=.GlobalEnv)
 on.exit({do.call(RNGkind,as.list(old_kind));if(had_seed)assign(".Random.seed",old_seed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
 RNGkind("L'Ecuyer-CMRG");set.seed(seed)
 seeds <- vector("list",simulations)
 seeds[[1]] <- get(".Random.seed",envir=.GlobalEnv)
 if(simulations>1L)for(b in 2:simulations)seeds[[b]] <- parallel::nextRNGStream(seeds[[b-1L]])
 workers <- min(workers,simulations)
 chunks <- split(seq_len(simulations),rep(seq_len(workers),length.out=simulations))
 nfeatures <- nrow(d);nsamples <- ncol(d)
 simulate_chunk <- function(indices) {
  exceed <- valid <- integer(nfeatures)
  for(i in seq_along(indices)) {
   b <- indices[i]
   assign(".Random.seed",seeds[[b]],envir=.GlobalEnv)
   simulated <- matrix(rnbinom(nfeatures*nsamples,mu=as.vector(null$fitted.values),size=rep(1/dispersion,nsamples)),nrow=nfeatures,dimnames=dimnames(d$counts))
   simulated_full <- glmFit(simulated,design,dispersion=dispersion,offset=offset,prior.count=0)
   simulated_null <- glmFit(simulated,null_design,dispersion=dispersion,offset=offset,prior.count=0)
   statistic <- pmax(0,simulated_null$deviance-simulated_full$deviance)
   good <- is.finite(statistic)
   valid[good] <- valid[good]+1L
   exceed[good] <- exceed[good]+as.integer(statistic[good]>=observed[good]-1e-10)
   if(indices[1]==1L&&(i==1L||i%%100L==0L||i==length(indices))) {
    cat("[MONTE CARLO]",basename(dirname(out)),"simulation",b,"of",simulations,"per feature\n");flush.console()
   }
  }
  list(exceed=exceed,valid=valid)
 }
 pieces <- if(workers==1L)list(simulate_chunk(chunks[[1]])) else parallel::mclapply(chunks,simulate_chunk,mc.cores=workers,mc.set.seed=FALSE)
 if(any(vapply(pieces,inherits,logical(1),"try-error")))stop("A Monte Carlo worker could not fit a simulated model")
 exceed <- Reduce(`+`,lapply(pieces,`[[`,"exceed"))
 valid <- Reduce(`+`,lapply(pieces,`[[`,"valid"))
 complete <- valid==simulations
 p <- (exceed+1)/(simulations+1)
 p[!complete] <- NA_real_
 lower <- ifelse(exceed==0,0,qbeta(.025,pmax(1,exceed),simulations-exceed+1))
 upper <- ifelse(exceed==simulations,1,qbeta(.975,exceed+1,pmax(1,simulations-exceed)))
 lower[!complete] <- upper[!complete] <- NA_real_
 data.frame(interaction_log2FC=drop(fit$coefficients %*% contrast)/log(2),LR=observed,
  PValue=p,MC_exceedances=exceed,MC_simulations=simulations,MC_valid_simulations=valid,
  MC_standard_error=sqrt(p*(1-p)/(simulations+1)),MC_tail_probability_CI95_lower=lower,
  MC_tail_probability_CI95_upper=upper,MC_minimum_pvalue=1/(simulations+1),MC_seed=seed,
  test_status=ifelse(complete,"tested","insufficient evidence: nonfinite simulated model statistics"),row.names=rownames(d))
}
