#' Person bootstrap for dynamic NIRA
#'
#' Resamples PERSONS with replacement, re-estimates J, h (and person effects)
#' with the original fitKIM() settings in every resample, and recomputes the
#' NIRA effects. The resulting intervals and rank frequencies therefore describe
#' sampling (estimation) uncertainty of the intervention targets, not Monte
#' Carlo variability; the latter is quantified separately by kimMCError().
#' The dose is held fixed at the original value (a fixed intervention); with
#' `recalibrate = TRUE` an equal-probability dose is recalibrated in every
#' resample instead. Selection and rank frequencies are descriptive, not p-values.
#' Percentile intervals for L1-penalised estimates can under-cover when effects
#' sit near zero; their coverage should be checked by simulation.
#' @param obj Output of runKIMNIRA.
#' @param data Original data.
#' @param n_boot Number of resamples. Default 2000. Small tail probabilities
#'   require many valid resamples; numerical resolution is not calibration.
#' @param seed Task seed.
#' @param parallel,ncores Use PSOCK parallel workers.
#' @param progress Print progress.
#' @param interval Percentile or bias-corrected intervals.
#' @param p_adjust Multiplicity adjustment for the reported p-values, passed to
#'   [stats::p.adjust] (default `"BH"`). Two families are adjusted separately:
#'   each node's effect against no change, and each node pair's benefit difference
#'   against equality. These p-values resample PERSONS, so unlike a test over
#'   Monte Carlo trajectories they do not shrink as the trajectory count grows.
#'   These approximate tests need coverage/type-I-error validation for the
#'   fitted procedure. They need not agree with unadjusted percentile/BC CIs.
#' @param p_method Centred bootstrap-error test (default), or an inclusive-tail
#'   percentile approximation. Degenerate nonzero centred tests return NA.
#' @param recalibrate Only for results run with `dose_scale = "probability"`:
#'   recalibrate the per-target doses in each resample (default `FALSE`, i.e. the
#'   original logit doses are held fixed).
#' @param checkpoint_dir Optional directory for restartable per-resample results.
#'   An input manifest prevents reuse with different data, settings or task seeds.
#' @export
kimStabilityBootstrap<-function(obj,data,n_boot=2000L,seed=2025L,
  parallel=FALSE,ncores=NULL,progress=TRUE,interval=c("percentile","bc"),
  p_adjust="BH",checkpoint_dir=NULL,p_method=c("centered","percentile"),
  recalibrate=FALSE){
  interval<-match.arg(interval)
  p_method<-match.arg(p_method)
  p_adjust<-match.arg(p_adjust,stats::p.adjust.methods)
  if(!inherits(obj,"kim_nira_result"))stop("obj must be a runKIMNIRA result.")
  fit<-obj$kim_fit;spec<-fit$settings
  if(is.null(spec))stop("Refit with current fitKIM to retain estimation settings.")
  if(length(n_boot)!=1L||!is.numeric(n_boot)||!is.finite(n_boot)||n_boot<2||
     n_boot>.Machine$integer.max||n_boot!=as.integer(n_boot))stop("n_boot must be an integer >= 2.")
  if(is.null(data[[spec$id]]))stop("Bootstrap data are missing the fitted person-ID column.")
  data<-data[as.character(data[[spec$id]]) %in% fit$person_ids,,drop=FALSE]
  if(!nrow(data))stop("No bootstrap rows match fitted person IDs.")
  set.seed(seed);task_seeds<-sample.int(.Machine$integer.max,n_boot)
  if(!is.null(checkpoint_dir)){
    dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
    checkpoint_dir<-normalizePath(checkpoint_dir,mustWork=TRUE)
    manifest<-list(version=1L,settings=spec,params=obj$params,delta=obj$delta,
      engine=obj$method,J=fit$J,h=fit$h,person_shifts=fit$person_shifts,
      person_J=fit$person_J,recalibrate=recalibrate,
      data=data,task_seeds=task_seeds)
    mp<-file.path(checkpoint_dir,'manifest.rds')
    if(file.exists(mp)){
      if(!identical(readRDS(mp),manifest))stop('Checkpoint inputs differ; use a new checkpoint directory.')
    }else saveRDS(manifest,mp)
  }
  one<-function(b){
    path<-if(is.null(checkpoint_dir))NULL else file.path(checkpoint_dir,sprintf('b%04d.rds',b))
    if(!is.null(path)&&file.exists(path)){
      cached<-tryCatch(readRDS(path),error=function(e)NULL)
      if(is.list(cached)&&length(cached$values)==length(fit$vars)+1L&&length(cached$error)==1L)return(cached)
      warning('Recomputing an unreadable or incomplete bootstrap checkpoint: ',basename(path))
    }
    set.seed(task_seeds[b]);ids<-unique(data[[spec$id]])
    sampled<-ids[sample.int(length(ids),length(ids),replace=TRUE)]
    boot<-do.call(rbind,lapply(seq_along(sampled),function(i){
      z<-data[data[[spec$id]]==sampled[i],,drop=FALSE];z[[spec$id]]<-paste0("B",i);z
    }))
    args<-spec;args$data<-boot;args$verbose<-FALSE;args$estimate_between<-FALSE
    args$seed<-as.integer(task_seeds[b] %% 1000000000L)
    # Per-resample warnings would fire hundreds of times; the flags are still
    # recorded on each refit and surfaced by the caller.
    args$slope_warn_singular<-FALSE
    result<-tryCatch({
      fb<-do.call(fitKIM,args)
      if(length(fb$node_fit_failed))stop("At least one node fit failed.")
      pa<-obj$params
      recal<-isTRUE(recalibrate)&&identical(pa$dose_scale,"probability")
      r<-runKIMNIRA(fb,intervention_type=pa$intervention_type,method=obj$method,
        delta=if(recal)NULL else obj$delta,
        dose_scale=if(recal)"probability" else "logit",target_dp=if(recal)pa$target_dp else NULL,
        calibrate_on=if(is.null(pa$calibrate_on))"stationary" else pa$calibrate_on,
        n_sim=if(obj$method=="mc")pa$n_sim else 1,n_steps=if(obj$method=="mc")pa$n_steps else 2,
        burn_in=if(obj$method=="mc")pa$burn_in else 0,
        paired=if(is.null(pa$paired)||is.na(pa$paired))TRUE else pa$paired,
        trajectory=FALSE,   # only effects/base_mean are used; paths would be 500x wasted memory
        diagnose=FALSE,     # multi-start pilot chains would be run 500x for nothing
        seed=as.integer(task_seeds[b] %% 1000000000L+1L),verbose=FALSE)
      list(values=c(r$effects,r$base_mean),error=NA_character_)
    },error=function(e)list(values=rep(NA_real_,length(fit$vars)+1),error=conditionMessage(e)))
    if(!is.null(path))saveRDS(result,path)
    result
  }
  if(parallel){
    if(is.null(ncores))ncores<-min(4L,parallel::detectCores())
    cl<-parallel::makePSOCKcluster(ncores);on.exit(parallel::stopCluster(cl),add=TRUE)
    env<-environment(runKIMNIRA)
    if(isNamespace(env)){
      # Every worker loads the installed package. A development session that
      # obtained the master namespace through pkgload::load_all() therefore has
      # to install the same version before running in parallel (or use
      # parallel = FALSE).
      ok <- tryCatch({parallel::clusterEvalQ(cl,library(KimNIRA));TRUE},
                     error=function(e)conditionMessage(e))
      if(!isTRUE(ok))
        stop("Parallel workers could not load an installed KimNIRA (",ok,
             "). Install the package, or use parallel = FALSE.",call.=FALSE)
      v <- unlist(parallel::clusterEvalQ(cl,as.character(utils::packageVersion("KimNIRA"))))
      master_v <- as.character(utils::packageVersion("KimNIRA"))
      if(any(v!=master_v))
        stop("Parallel workers loaded a different KimNIRA version (",
             paste(unique(v),collapse=", "),") than the master (",master_v,
             "). Install the same version, or use parallel = FALSE.",call.=FALSE)
    }else{
      nms<-ls(env,all.names=TRUE)
      nms<-nms[vapply(nms,function(nm)is.function(get(nm,env)),logical(1))]
      parallel::clusterExport(cl,nms,envir=env)
      parallel::clusterEvalQ(cl,.kim_state_cache<-new.env(parent=emptyenv()))
    }
    ans<-parallel::parLapply(cl,seq_len(n_boot),one)
  }else{
    ans<-lapply(seq_len(n_boot),function(b){
      z<-one(b);if(progress&&(b%%25==0||b==n_boot))message("Bootstrap ",b,"/",n_boot);z
    })
  }
  errors<-vapply(ans,`[[`,character(1),"error");all<-do.call(rbind,lapply(ans,`[[`,"values"))
  ok<-complete.cases(all)
  if(sum(ok)<2){
    first <- unique(errors[!ok]); first <- first[seq_len(min(5L,length(first)))]
    stop("Fewer than two valid bootstrap resamples. First errors:
",
         paste(first,collapse="
"),call.=FALSE)
  }
  if(any(!ok))warning(sum(!ok)," bootstrap resamples failed; inspect failures.")
  eff<-all[ok,seq_along(fit$vars),drop=FALSE];colnames(eff)<-fit$vars
  base<-all[ok,ncol(all)];p<-ncol(eff)
  ci<-function(x,point)if(interval=="bc").kim_bc_int(x,point)else stats::quantile(x,c(.025,.975),names=FALSE)
  ints<-t(vapply(seq_len(p),function(j)ci(eff[,j],obj$effects[j]),numeric(2)))
  benefit<-if(obj$params$intervention_type=="alleviating")-eff else eff
  rankprops<-matrix(0,p,p,dimnames=list(fit$vars,paste0("rank_",seq_len(p))))
  for(r in seq_len(nrow(eff))){
    values<-sort(unique(benefit[r,]),decreasing=TRUE);pos<-1L
    for(v in values){
      nodes<-which(benefit[r,]==v);m<-length(nodes)
      rankprops[nodes,pos:(pos+m-1L)]<-rankprops[nodes,pos:(pos+m-1L)]+1/m;pos<-pos+m
    }
  }
  rankprops<-rankprops/nrow(eff)
  st<-data.frame(node=fit$vars,effect=unname(obj$effects),se_effect=apply(eff,2,stats::sd),
    ci_effect_lower=ints[,1],ci_effect_upper=ints[,2],percenttop_1=rankprops[,1])
  # Significance is assessed over the PERSON bootstrap, so it does not shrink as the
  # Monte Carlo trajectory count grows. A test whose sample size is n_sim asks only
  # whether the simulated difference is exactly zero and marks every node `***`.
  p_zero<-vapply(seq_len(p),function(j).kim_boot_p(eff[,j],
    point=obj$effects[j],method=p_method),numeric(1))
  st$p_value<-p_zero
  st$p_adjusted<-stats::p.adjust(p_zero,method=p_adjust)
  if(anyNA(p_zero))warning("Some bootstrap distributions are degenerate; their p-values are NA and no stars are drawn.",call.=FALSE)
  # The finite-resample p-value floor depends on the selected test. Warn when it
  # what denies the three-star level: at the default n_boot it is not, and warning
  # on every run would be noise rather than information.
  p_floor<-if(p_method=="centered")1/(sum(ok)+1) else 2/(sum(ok)+1)
  if(p_floor>0.001&&any(is.finite(p_zero)&p_zero<=p_floor*1.0000001))
    warning("n_boot = ",n_boot," cannot produce a p-value below ",signif(p_floor,3),
      ", so some nodes are capped at '",.kim_star_label(p_floor),
      "'. More resamples improve tail resolution, not statistical validity.",call.=FALSE)
  post<-sweep(eff,1,base,"+")
  post_ci<-t(vapply(seq_len(p),function(j)ci(post[,j],obj$base_mean+obj$effects[j]),numeric(2)))
  bci<-ci(base,obj$base_mean)
  plot_data<-data.frame(node=c("baseline",fit$vars),mean=c(obj$base_mean,obj$base_mean+obj$effects),
    ciLower=c(bci[1],post_ci[,1]),ciUpper=c(bci[2],post_ci[,2]))
  pairwise<-data.frame()
  if(p>1){
    ij<-utils::combn(p,2)
    diffs<-vapply(seq_len(ncol(ij)),function(i)benefit[,ij[1,i]]-benefit[,ij[2,i]],numeric(nrow(eff)))
    if(is.null(dim(diffs)))diffs<-matrix(diffs,ncol=1)
    g<-if(obj$params$intervention_type=="alleviating")-obj$effects else obj$effects
    point<-g[ij[1,]]-g[ij[2,]]
    critical<-unname(stats::quantile(apply(abs(sweep(diffs,2,point,"-")),1,max),.95))
    pairci<-apply(diffs,2,stats::quantile,c(.025,.975))
    p_pair<-vapply(seq_len(ncol(diffs)),function(j).kim_boot_p(diffs[,j],
      point=point[j],method=p_method),numeric(1))
    pairwise<-data.frame(node1=fit$vars[ij[1,]],node2=fit$vars[ij[2,]],benefit_difference=point,
      lower=pairci[1,],upper=pairci[2,],simultaneous_lower=point-critical,simultaneous_upper=point+critical,
      p_value=p_pair,p_adjusted=stats::p.adjust(p_pair,method=p_adjust))
  }
  out<-list(effects=obj$effects,base_mean=obj$base_mean,base_ci_lower=bci[1],base_ci_upper=bci[2],
    resample_base=base,resample_effects=eff,stat=st,pairwise=pairwise,plot_data=plot_data,
    rank_proportions=rankprops,stable_ranking=obj$ranked_targets,
    kim_fit=fit,source_signature=.kim_inference_signature(obj),
    failures=data.frame(replicate=which(!ok),message=errors[!ok]),
    parameters=list(engine=obj$method,intervention_type=obj$params$intervention_type,
      delta=obj$delta,burn_in=obj$params$burn_in,interval=interval,p_adjust=p_adjust,
      p_method=p_method,inference_version=2L,p_floor=p_floor,
      n_boot=n_boot,n_valid=sum(ok),seed=seed))
  class(out)<-"kim_stability_bootstrap";out
}
#' @export
print.kim_stability_bootstrap<-function(x,...){
  cat("Person bootstrap: ",x$parameters$n_valid,"/",x$parameters$n_boot," valid resamples\n",sep="")
  print(x$stat,row.names=FALSE);invisible(x)
}
