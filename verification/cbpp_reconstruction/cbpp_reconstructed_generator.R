# Independent candidate reconstruction of the upstream CBPP grouped-binomial RR calibration.
# NOT recovered historical generator. Original archived code is absent.
# Run from package root; uses fixed per-replicate seeds, checkpointed batches, 6 processes.
options(warn=1)
stopifnot(requireNamespace('lme4',quietly=TRUE),
          requireNamespace('statmod',quietly=TRUE),
          requireNamespace('dplyr',quietly=TRUE))
suppressPackageStartupMessages({library(lme4);library(statmod);library(dplyr)})
source('code/functions_05_empirical_additional_examples.R',local=TRUE)
mode <- tolower(Sys.getenv('RUN_MODE','full'))
stopifnot(mode %in% c('full','pilot'))
B <- if (mode=='full') 1000L else 12L
workers <- as.integer(Sys.getenv('WORKERS','6'))
if(is.na(workers)||workers<1L) workers <- 1L
workers <- min(workers,6L)
output <- file.path('results','cbpp_reconstruction',mode)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
base <- read.csv('reference/cbpp_RR_observed_residuals.csv',check.names=FALSE)
stopifnot(nrow(base)==56L,all(c('herd','incidence','size','period','period_num')%in% names(base)))
base$herd <- factor(base$herd)
base$period <- factor(base$period)
base$period_num <- as.numeric(base$period_num)
base <- base[order(base$herd,base$period_num),]
base$obs <- seq_len(nrow(base))
formula <- cbind(incidence,size-incidence) ~ period+(1|herd)
ctrl <- glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e5))
fit <- glmer(formula,data=base,family=binomial(),control=ctrl,nAGQ=1)
refdh <- read.csv('reference/cbpp_DHARMa_AIC_B1000.csv')
obs_aic <- as.numeric(AIC(fit))
ref_aic <- refdh$AIC_baseline[1]
cat(sprintf('CBPP_BASELINE_AIC %.12f ARCHIVED %.12f\n',obs_aic,ref_aic))
if(abs(obs_aic-ref_aic)>1e-4) stop('Baseline AIC does not match; do not run reconstructed calibration')
fit_olre <- glmer(cbind(incidence,size-incidence)~period+(1|herd)+(1|obs),
                  data=base,family=binomial(),control=ctrl,nAGQ=1)
olre_aic <- AIC(fit_olre)
cat(sprintf('CBPP_OLRE_AIC %.12f ARCHIVED %.12f\n',olre_aic,refdh$AIC_OLRE[1]))
if(abs(olre_aic-refdh$AIC_OLRE[1])>1e-3) stop('OLRE AIC does not match')
# Sequential binomial predictive PIT integrating each herd's Gaussian RI by quadrature.
# This is a transparent independently reconstructed candidate; it may differ from historical method.
rr_cbpp <- function(fitted, dat, seed) {
  set.seed(seed)
  dat <- dat[order(dat$herd,dat$period_num),,drop=FALSE]
  b <- lme4::fixef(fitted)
  X <- model.matrix(~period,data=dat)
  stopifnot(identical(colnames(X),names(b)))
  sig <- as.numeric(attr(lme4::VarCorr(fitted)$herd,'stddev'))[1]
  quad <- statmod::gauss.quad.prob(50,dist='normal')
  nodes <- quad$nodes*sig
  prob0 <- as.numeric(quad$weights)
  res <- rep(NA_real_,nrow(dat)); U <- res; pred <- res
  for (id in levels(dat$herd)) {
    rr <- which(dat$herd==id)
    lw <- rep(0,length(nodes))
    for(j in rr) {
      w <- exp(log(prob0)+lw-max(log(prob0)+lw))
      w <- w/sum(w)
      p <- plogis(drop(X[j,,drop=FALSE]%*%b)+nodes)
      y <- dat$incidence[j]; m <- dat$size[j]
      cdf_lo <- sum(w*pbinom(y-1L,size=m,prob=p))
      cdf_hi <- sum(w*pbinom(y,size=m,prob=p))
      U[j] <- cdf_lo+runif(1)*(cdf_hi-cdf_lo)
      res[j] <- qnorm(pmin(1-.Machine$double.eps,pmax(.Machine$double.eps,U[j])))
      pred[j] <- sum(w*p)
      lw <- lw+dbinom(y,size=m,prob=p,log=TRUE)
    }
  }
  dat$rr <- res;dat$u <- U;dat$pred_prob <- pred;dat$pred_mean <- dat$size*pred
  dat
}
index_vars <- c('period_num','pred_prob','size')
obsrr <- rr_cbpp(fit,base,2026L)
obsstats <- component_stats(obsrr,index_vars,'herd','period_num')
write.csv(obsrr,file.path(output,'cbpp_RR_observed_residuals.csv'),row.names=FALSE)
write.csv(obsstats,file.path(output,'cbpp_observed_stats.csv'),row.names=FALSE)
refobs <- read.csv('reference/cbpp_RR_B1000_component_calibration.csv')
refobs <- refobs[,c('statistic','observed')]
obs_diff <- vapply(refobs$statistic,function(nm) {
 if(!(nm %in% names(obsstats))) return(Inf)
 abs(obsstats[[nm]][1]-refobs$observed[refobs$statistic==nm][1])
},numeric(1))
write.csv(data.frame(statistic=refobs$statistic,absolute_difference=obs_diff),
          file.path(output,'observed_stats_differences.csv'),row.names=FALSE)
cat('CBPP_OBSERVED_MAX_ABS_DIFF',max(obs_diff),'\n')
if(mode == 'full' && (!all(is.finite(obs_diff)) || max(obs_diff) > 1e-7)) {
  stop('Observed RR statistics fail archived reference at 1e-7; full 1000-refit run halted to avoid unjustified reproduction claim. Run pilot for exploratory checks.')
}
# Explicit gate: do NOT claim equality if different.
sim_one <- function(i) {
  tryCatch({
    y <- simulate(fit,nsim=1,seed=10000L+i)[[1]]
    if(is.matrix(y)) y <- y[,1]
    y <- as.numeric(y)
    if(length(y)!=nrow(base)||any(y<0)||any(y>base$size)||any(abs(y-round(y))>1e-8)) stop('Unexpected grouped binomial simulated response')
    dat <- base;dat$incidence <- y
    f <- suppressWarnings(glmer(formula,data=dat,family=binomial(),control=ctrl,nAGQ=1))
    r <- rr_cbpp(f,dat,20000L+i)
    s <- component_stats(r,index_vars,'herd','period_num')
    s$replicate <- i
    s$singular <- lme4::isSingular(f)
    s$failure <- FALSE
    s
  },error=function(e) data.frame(replicate=i,singular=NA, failure=TRUE,error=conditionMessage(e)))
}
progress <- file.path(output,'checkpoint.csv')
if(file.exists(progress)) {
 prev <- read.csv(progress,check.names=FALSE)
 prev <- prev[!duplicated(prev$replicate),,drop=FALSE]
 if(!all(c(names(obsstats),'replicate','singular','failure')%in%names(prev))) stop('Existing checkpoint schema incompatible')
} else prev <- data.frame()
done <- if(nrow(prev)) prev$replicate else integer(0)
remaining <- setdiff(seq_len(B),done)
cat('CBPP_RUN',mode,'replicates',B,'already_completed',length(done),'workers',workers,'\n')
if (length(remaining)>0L) for (block_start in seq(1L,length(remaining),by=24L)) {
 idx <- remaining[seq.int(block_start,min(block_start+23L,length(remaining)))]
 if(length(idx)==0) next
 if(.Platform$OS.type=='unix'&&workers>1) {
   rows <- parallel::mclapply(idx,sim_one,mc.cores=workers,mc.set.seed=FALSE)
 } else rows <- lapply(idx,sim_one)
 new <- dplyr::bind_rows(rows)
 prev <- dplyr::bind_rows(prev,new)
 prev <- prev[order(prev$replicate),,drop=FALSE]
 tmp <- paste0(progress,'.tmp')
 write.csv(prev,tmp,row.names=FALSE)
 if(!file.rename(tmp,progress)) stop('Could not atomically update checkpoint')
 cat('CBPP_CHECKPOINT',nrow(prev),'/',B,'failed',sum(prev$failure,na.rm=TRUE),'\n')
}
if(nrow(prev)!=B) stop('Incomplete CBPP checkpoint')
if(any(prev$failure,na.rm=TRUE)) {
 cat('CBPP_REFIT_FAILURES',sum(prev$failure,na.rm=TRUE),'see checkpoint; no exact certification\n')
}
valid <- prev[!prev$failure, ,drop=FALSE]
simstats <- valid[,names(obsstats),drop=FALSE]
write.csv(simstats,file.path(output,'cbpp_RR_B1000_simulated_stats.csv'),row.names=FALSE)
if(nrow(simstats)>=50) {
 cal <- calibrate_components(obsstats,simstats,0.05)
 write.csv(cal,file.path(output,'cbpp_RR_component_calibration.csv'),row.names=FALSE)
 cat('CBPP_RECONSTRUCTED_P_COMBO',cal$p_combo_minP[1],'ARCHIVED',read.csv('reference/cbpp_RR_B1000_component_calibration.csv')$p_combo_minP[1],'\n')
}
write.csv(data.frame(B_requested=B,B_success=nrow(valid),failures=sum(prev$failure,na.rm=TRUE),singular=sum(valid$singular,na.rm=TRUE),
                     AIC_baseline=obs_aic,AIC_OLRE=olre_aic,max_observed_stat_abs_diff=max(obs_diff)),
          file.path(output,'run_summary.csv'),row.names=FALSE)
cat('CBPP_CANDIDATE_GENERATOR_RUN_COMPLETE; NOT HISTORICAL GENERATOR CERTIFICATION\n')
