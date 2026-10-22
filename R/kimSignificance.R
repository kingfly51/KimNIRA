#' Intervention effects with bootstrap intervals
#' @param boot Output of kimStabilityBootstrap.
#' @param stars Significance labels from the person-bootstrap p-values: `"vs_zero"`
#'   for an effect distinguishable from no change, `"vs_best"` appending `"="` when
#'   a node cannot be distinguished from the top-ranked target, `"both"` for both.
#'   `"auto"` uses `"vs_zero"` when verified p-values are present.
#'   These approximate tests need not agree with unadjusted effect intervals.
#' @return A ggplot object showing point estimates and percentile or BC intervals.
#' @export
kimDeltaPlot<-function(boot,stars=c("auto","both","none","vs_zero","vs_best")){
 if(!inherits(boot,"kim_stability_bootstrap"))
  stop("`boot` must be the output of kimStabilityBootstrap(). ",
       "For runKIMNIRA output use kimPlotNIRA() instead.")
 stars<-match.arg(stars)
 if(stars=="auto")stars<-if("p_adjusted"%in%names(boot$stat))"vs_zero" else "none"
 st<-boot$stat
 st<-st[order(st$effect,decreasing=boot$parameters$intervention_type=="aggravating"),]
 st$node<-factor(st$node,levels=st$node)
 st$star<-if(stars=="none")rep("",nrow(st)) else .kim_star_text(boot,as.character(st$node),stars)
 st$star_y<-ifelse(st$effect>=0,st$ci_effect_upper,st$ci_effect_lower)
 labelled<-!is.na(st$star)&nzchar(st$star)
 gg<-ggplot2::ggplot(st,ggplot2::aes(x=.data$node,y=.data$effect))+
  ggplot2::geom_hline(yintercept=0,linetype=2,colour="grey40")+
  ggplot2::geom_col(ggplot2::aes(fill=.data$effect >= 0),width=.65,alpha=.85)+
  ggplot2::scale_fill_manual(values=c("FALSE"="#1B9E77","TRUE"="#D95F02"),guide="none")+
  ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$ci_effect_lower,ymax=.data$ci_effect_upper),width=.2)
 if(any(labelled)){
  # Two layers with scalar vjust: vjust is not an aesthetic in older ggplot2.
  up<-labelled&st$effect>=0;dn<-labelled&st$effect<0
  if(any(up))gg<-gg+ggplot2::geom_text(data=st[up,,drop=FALSE],
   ggplot2::aes(x=.data$node,y=.data$star_y,label=.data$star),vjust=-0.7,size=4)
  if(any(dn))gg<-gg+ggplot2::geom_text(data=st[dn,,drop=FALSE],
   ggplot2::aes(x=.data$node,y=.data$star_y,label=.data$star),vjust=1.7,size=4)
 }
 note<-if(stars=="vs_best")
  "= not distinguishable from the top target (person bootstrap, adjusted)" else
  if(stars=="vs_zero")
  "* p<.05, ** p<.01, *** p<.001 for the effect vs no change (person bootstrap, adjusted)" else
  paste0("* p<.05, ** p<.01, *** p<.001 for the effect vs no change; ",
         " = not distinguishable from the top target (person bootstrap, adjusted)")
 gg+ggplot2::scale_y_continuous(labels=function(x) format(100*x,trim=TRUE))+
  ggplot2::labs(x=NULL,y="Change in mean activation (percentage points)",
   title="Intervention effect with 95% bootstrap CI",
   caption=if(any(labelled))paste0(note,"; method: ",boot$parameters$p_method,
     "; adjustment: ",boot$parameters$p_adjust,". CI and test need not agree.") else NULL)+
  ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
}
