# Significance labels for a set of nodes, from a bootstrap result's person-level
# p-values. `vs_zero` marks an effect distinguishable from no change; `vs_best`
# appends "=" when a node is NOT distinguishable from the top-ranked target, which
# is the decision-relevant pairwise statement. Returns "" where nothing applies.
.kim_star_text <- function(boot, nodes, stars = c("both", "none", "vs_zero", "vs_best")) {
  stars <- match.arg(stars)
  out   <- rep("", length(nodes))
  if (stars == "none" || is.null(boot$stat) || !"p_adjusted" %in% names(boot$stat))
    return(out)
  if (!identical(boot$parameters$inference_version, 2L))
    stop("Legacy/unverified bootstrap p-values: rerun kimStabilityBootstrap before drawing stars.")
  if (stars %in% c("vs_zero", "both"))
    out <- .kim_star_label(boot$stat$p_adjusted[match(nodes, boot$stat$node)])
  if (stars %in% c("vs_best", "both")) {
    pw   <- boot$pairwise
    best <- boot$stable_ranking[1]
    if (!is.null(pw) && nrow(pw) > 0L && "p_adjusted" %in% names(pw) && !is.null(best)) {
      for (i in seq_along(nodes)) {
        if (is.na(nodes[i]) || identical(nodes[i], best)) next
        row <- which((pw$node1 == best & pw$node2 == nodes[i]) |
                     (pw$node2 == best & pw$node1 == nodes[i]))
        if (length(row) && isTRUE(pw$p_adjusted[row[1]] > .05))
          out[i] <- paste0(out[i], "=")
      }
    }
  }
  out
}

# Match inference to a specific fitted model, dose and numerical run.
.kim_inference_signature <- function(obj) {
  f <- obj$kim_fit
  list(vars=f$vars,J=f$J,h=f$h,person_shifts=f$person_shifts,person_J=f$person_J,
       person_ids=f$person_ids,settings=f$settings,method=obj$method,
       params=obj$params,delta=obj$delta,effects=obj$effects,base_mean=obj$base_mean)
}

# Preserve original colours, extend short palettes, and respect named overrides.
.kim_plot_colors <- function(ids, colors = NULL, defaults) {
  n <- length(ids)
  out <- stats::setNames(rep(NA_character_, n), ids)
  if (length(colors)) {
    if (anyNA(colors)) stop("`colors` cannot contain NA.")
    if (!is.null(names(colors))) {
      if (any(!nzchar(names(colors))) || anyDuplicated(names(colors)) ||
          any(!names(colors) %in% ids)) stop("Named colors must have unique, valid node/rank names.")
      out[names(colors)] <- colors
    } else out[seq_len(min(n, length(colors)))] <- colors[seq_len(min(n, length(colors)))]
  }
  # Spread additional hues around the wheel rather than taking adjacent hues.
  extra <- grDevices::hcl.colors(max(30L, 2L*n), "Dark 3")
  idx <- unique(floor((((seq_along(extra)-1)*0.61803398875) %% 1)*length(extra))+1L)
  pool <- unique(c(defaults, extra[c(idx, setdiff(seq_along(extra),idx))]))
  available <- pool[!pool %in% out[!is.na(out)]]
  out[is.na(out)] <- available[seq_len(sum(is.na(out)))]
  grDevices::col2rgb(out)
  out
}
