#' Deprecated functions from ICONS 0.1.x
#'
#' These keep scripts written against ICONS 0.1.x running.  Each returns the
#' old shape of output and warns once per session.  They will be removed in a
#' future release.
#'
#' | 0.1.x | 0.2.0 |
#' | --- | --- |
#' | `dense()` | [icons_detect()] |
#' | `param_tuning_sigmau()` | [icons_tune()] |
#' | `k.elbow()` | [n_factors()] |
#' | `plotMatrix()` | [plot_matrix()] |
#' | `get_membership()` | [as_membership()] |
#' | `get_index()` | [block_index()] |
#' | `get_vectorform()` | [half_vec()] |
#'
#' `scfa()` kept its name but changed signature: it now takes a partition
#' object (or membership vector) instead of the `CID`/`Clist` pair, returns
#' the closed-form maximum likelihood estimators of Yang, Ma, Bi and Chen
#' (2024) rather than `cov(F_hat)`, and no longer builds `p` by `p` residual
#' matrices.  See `vignette("ICONS")` for the migration.
#'
#' @param W_original,Wp Adjacency matrix.
#' @param threshold,lambda Detection settings.
#' @param data Data matrix.
#' @param CID Community sizes.
#' @param Clist Reordered variable indices.
#' @param prctile_vec Percentiles (0-100) for the threshold grid.
#' @param lam_vec Lambda grid.
#' @param method Ignored; retained for signature compatibility.
#' @param ncores,use_parallel Parallel settings.
#' @param k Maximum number of factors.
#' @param epsilon Ignored; the closed-form estimators need no regularisation.
#' @param z Matrix to plot.
#' @param blockid Block numbers.
#' @param dist Symmetric matrix.
#' @param filepath,save.image,width,height,res,format,colorbar.range Plot
#'   arguments.
#' @param ... Passed through.
#'
#' @return The same shapes ICONS 0.1.x returned.
#'
#' @name ICONS-deprecated
#' @keywords internal
NULL

warn_once <- local({
  seen <- new.env(parent = emptyenv())
  function(old, new) {
    if (is.null(seen[[old]])) {
      assign(old, TRUE, envir = seen)
      warning(old, "() is deprecated in ICONS 0.2.0; use ", new, "() instead.",
              call. = FALSE)
    }
  }
})

# Rebuild a 0.1.x-style partition (CID / Clist) as an icons_partition.  The old
# convention appended every unassigned variable as a trailing pseudo-community,
# which is why `scfa()` had a `remove.singletons` argument.
legacy_partition <- function(CID, Clist, p) {
  CID <- as.integer(CID)
  Clist <- as.integer(Clist)
  # 0.1.x sized this vector by length(Clist) alone, so an index beyond that
  # left NA holes; size it to cover every index actually referenced.
  membership <- integer(max(p, Clist, 0L))
  end <- cumsum(CID)
  start <- end - CID + 1L
  for (k in seq_along(CID)) {
    membership[Clist[start[k]:end[k]]] <- k
  }
  membership
}

#' @rdname ICONS-deprecated
#' @export
dense <- function(W_original, threshold = 0.5, lambda = 0.5) {
  warn_once("dense", "icons_detect")
  part <- icons_detect(W_original, threshold = threshold, lambda = lambda)
  W <- as_weight_matrix(W_original, arg = "W_original")
  W <- W - diag(diag(W))
  list(W_dense = W[part$order, part$order],
       Clist = part$order,
       CID = c(part$sizes, if (part$n_singletons > 0L) part$n_singletons))
}

#' @rdname ICONS-deprecated
#' @export
param_tuning_sigmau <- function(Wp, data, prctile_vec, lam_vec,
                                method = "Sample", ncores = NULL,
                                use_parallel = TRUE) {
  warn_once("param_tuning_sigmau", "icons_tune")
  if (is.null(ncores)) ncores <- 1L
  if (!use_parallel) ncores <- 1L
  res <- icons_tune(Wp, data, probs = prctile_vec / 100, lambda = lam_vec,
                    criterion = "legacy", ncores = ncores)
  all <- data.frame(Lambda = res$grid$lambda, CutOff = res$grid$threshold,
                    SigmaU = res$grid$criterion)
  list(lambda_out = res$lambda, cut_out = res$threshold, all = all)
}

#' @rdname ICONS-deprecated
#' @export
k.elbow <- function(data, CID, Clist, k = length(CID), method = "Sample",
                    epsilon = 1e-6) {
  warn_once("k.elbow", "n_factors")
  p <- ncol(as_data_matrix(data, arg = "data"))
  memb <- legacy_partition(CID, Clist, p)
  nf <- n_factors(data, memb, k_max = min(k, max(memb)), criterion = "legacy")
  nf$path$criterion
}

#' @rdname ICONS-deprecated
#' @export
plotMatrix <- function(z, filepath = NULL, save.image = FALSE,
                       width = 2500, height = 2350, res = 300,
                       format = "tiff", colorbar.range = NULL, ...) {
  warn_once("plotMatrix", "plot_matrix")
  plot_matrix(z, zlim = colorbar.range,
              file = if (isTRUE(save.image)) {
                if (is.null(filepath)) file.path(getwd(), "plot.tiff") else filepath
              } else NULL,
              format = format, width = width, height = height, res = res)
}

#' @rdname ICONS-deprecated
#' @export
get_membership <- function(CID, Clist) {
  warn_once("get_membership", "as_membership")
  legacy_partition(CID, Clist, length(Clist))
}

#' @rdname ICONS-deprecated
#' @export
get_index <- function(blockid, CID, Clist) {
  warn_once("get_index", "block_index")
  blockid <- as.integer(blockid)
  if (any(blockid < 1L | blockid > length(CID))) {
    abort_icons("One or more block IDs are out of range.",
                class = "icons_value_error")
  }
  end <- cumsum(CID)
  start <- end - CID + 1L
  list(indices = unlist(lapply(blockid, function(j) Clist[start[j]:end[j]]),
                        use.names = FALSE),
       block_ids = rep.int(blockid, CID[blockid]))
}

#' @rdname ICONS-deprecated
#' @export
get_vectorform <- function(dist) {
  warn_once("get_vectorform", "half_vec")
  half_vec(dist)
}
