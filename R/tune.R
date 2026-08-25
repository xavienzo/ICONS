#' Tune the detection threshold and shrinkage exponent
#'
#' Grid search over the edge threshold and \eqn{\lambda}.  Every candidate pair
#' is run through [icons_detect()] and [scfa()], and the pair minimising the
#' fit criterion is returned.
#'
#' The threshold matters far more than \eqn{\lambda}: it decides which edges
#' the peeling ever sees, so it is specified as a grid of quantiles of the
#' observed edge weights, which adapts to the scale of the data.
#'
#' @section Interpreting the result:
#' The criterion is evaluated on the same matrix used to select the parameters,
#' so it measures fit, not out-of-sample performance, and the minimised value
#' is optimistically biased.  Treat the selected `threshold`/`lambda` as a
#' data-driven default and the criterion surface --- `plot()` on the result ---
#' as the more informative output: a flat surface means the choice is not
#' well identified.
#'
#' @param W A square symmetric numeric matrix, typically `cor(data)`.
#' @param data The numeric data matrix `W` was computed from, observations in
#'   rows.
#' @param probs Numeric vector of quantiles of the edge weights to use as
#'   thresholds.  Defaults to `seq(0.90, 0.99, by = 0.01)`.
#' @param lambda Numeric vector of shrinkage exponents.  Defaults to
#'   `seq(0.4, 0.8, length.out = 5)`.
#' @param criterion Fit criterion, `"frobenius"` (default) or `"legacy"`; see
#'   [scfa_criterion()].
#' @param min_size Smallest retained community, passed to [icons_detect()].
#' @param ncores Number of worker processes.  Defaults to
#'   `getOption("ICONS.ncores", 1L)`; `1` runs sequentially.
#' @param verbose Logical; report progress.
#'
#' @return An object of class `"icons_tune"`:
#'   \describe{
#'     \item{`best`}{One-row data frame with the selected settings.}
#'     \item{`grid`}{Data frame with one row per candidate: `probs`,
#'       `threshold`, `lambda`, `k` (communities found), `n_singletons`,
#'       `criterion`, `relative`.}
#'     \item{`threshold`, `lambda`}{The selected values, for passing straight
#'       to [icons_detect()].}
#'   }
#'
#' @seealso [icons_detect()], [scfa()], [n_factors()]
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#' tuned <- icons_tune(W, sim, probs = c(0.90, 0.95), lambda = c(0.5, 0.7))
#' tuned
#'
#' part <- icons_detect(W, threshold = tuned$threshold, lambda = tuned$lambda)
#'
#' @export
icons_tune <- function(W,
                       data,
                       probs = seq(0.90, 0.99, by = 0.01),
                       lambda = seq(0.4, 0.8, length.out = 5),
                       criterion = c("frobenius", "legacy"),
                       min_size = 2L,
                       ncores = getOption("ICONS.ncores", 1L),
                       verbose = FALSE) {
  cl <- match.call()
  criterion <- match.arg(criterion)
  W <- as_weight_matrix(W, arg = "W")
  X <- as_data_matrix(data, arg = "data")
  if (ncol(X) != ncol(W)) {
    abort_icons("`data` has ", ncol(X), " columns but `W` is ", ncol(W),
                " by ", ncol(W), ".", class = "icons_dim_error")
  }
  if (!is.numeric(probs) || !length(probs) || any(probs <= 0 | probs >= 1)) {
    abort_icons("`probs` must be numbers strictly between 0 and 1.",
                class = "icons_value_error")
  }
  if (!is.numeric(lambda) || !length(lambda) || any(lambda < 0 | lambda > 1)) {
    abort_icons("`lambda` must be numbers between 0 and 1.",
                class = "icons_value_error")
  }
  min_size <- check_count(min_size, "min_size", lower = 2L)
  ncores <- check_count(ncores, "ncores", lower = 1L)

  thresholds <- unname(stats::quantile(W[upper.tri(W)], sort(probs)))
  grid <- expand.grid(lambda = sort(lambda), i = seq_along(thresholds),
                      KEEP.OUT.ATTRS = FALSE)
  grid$probs <- sort(probs)[grid$i]
  grid$threshold <- thresholds[grid$i]
  grid$i <- NULL

  one <- function(r) {
    part <- icons_detect(W, threshold = grid$threshold[r],
                         lambda = grid$lambda[r], min_size = min_size)
    if (!length(part$sizes)) {
      return(c(k = 0, n_singletons = ncol(W), criterion = NA_real_,
               relative = NA_real_))
    }
    fit <- tryCatch(
      suppressWarnings(scfa(X, part, criterion = criterion)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(k = length(part$sizes), n_singletons = part$n_singletons,
               criterion = NA_real_, relative = NA_real_))
    }
    c(k = fit$K, n_singletons = fit$n_singletons,
      criterion = fit$criterion$value, relative = fit$criterion$relative)
  }

  idx <- seq_len(nrow(grid))
  if (verbose) {
    message("ICONS: evaluating ", length(idx), " parameter combinations on ",
            ncores, " core(s)")
  }

  res <- if (ncores > 1L) {
    run_parallel(idx, one, ncores)
  } else {
    lapply(idx, one)
  }

  out <- cbind(grid[, c("probs", "threshold", "lambda")],
               as.data.frame(do.call(rbind, res)))
  out <- out[order(out$probs, out$lambda), ]
  rownames(out) <- NULL

  if (all(is.na(out$criterion))) {
    abort_icons("no parameter combination produced a fittable model; ",
                "try lower `probs` or a smaller `min_size`.",
                class = "icons_value_error")
  }
  best <- out[which.min(out$criterion), , drop = FALSE]

  structure(
    list(best = best, grid = out,
         threshold = best$threshold, lambda = best$lambda,
         probs = best$probs, criterion = criterion, call = cl),
    class = "icons_tune"
  )
}

# Cross-platform worker pool.  forking is used where it is available because it
# avoids re-serialising W, which can be large.
run_parallel <- function(idx, fun, ncores) {
  ncores <- min(ncores, length(idx))
  if (.Platform$OS.type != "windows") {
    res <- parallel::mclapply(idx, fun, mc.cores = ncores)
    failed <- vapply(res, inherits, logical(1L), what = "try-error")
    if (any(failed)) {
      stop("parallel evaluation failed: ",
           conditionMessage(attr(res[[which(failed)[1L]]], "condition")),
           call. = FALSE)
    }
    return(res)
  }
  cl <- parallel::makePSOCKcluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, requireNamespace("ICONS", quietly = TRUE))
  parallel::parLapplyLB(cl, idx, fun)
}

#' @export
print.icons_tune <- function(x, digits = 4L, ...) {
  cat("<ICONS parameter tuning>\n")
  cat("  grid       :", nrow(x$grid), "combinations\n")
  cat("  criterion  :", x$criterion, "\n")
  cat("  selected   : threshold =", format(x$threshold, digits = digits),
      "(quantile", format(x$probs, digits = 3), "), lambda =",
      format(x$lambda, digits = digits), "\n")
  cat("  gives      :", x$best$k, "communities,", x$best$n_singletons,
      "singletons, criterion", format(x$best$criterion, digits = digits), "\n")
  invisible(x)
}

#' @export
summary.icons_tune <- function(object, ...) {
  print(object)
  cat("\n")
  print(object$grid, row.names = FALSE, digits = 4L)
  invisible(object)
}

#' Plot the tuning criterion surface
#'
#' One line per `lambda`, criterion against the threshold quantile.  A flat
#' bundle of lines means the criterion does not really distinguish the
#' candidates.
#'
#' @param x An `"icons_tune"` object.
#' @param ... Passed to [graphics::matplot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' data(sim)
#' tuned <- icons_tune(cor(sim), sim, probs = c(0.9, 0.95), lambda = c(0.5, 0.7))
#' plot(tuned)
#'
#' @export
plot.icons_tune <- function(x, ...) {
  g <- x$grid
  lam <- sort(unique(g$lambda))
  pr <- sort(unique(g$probs))
  Z <- matrix(NA_real_, length(pr), length(lam))
  Z[cbind(match(g$probs, pr), match(g$lambda, lam))] <- g$criterion

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  graphics::matplot(pr, Z, type = "b", pch = 16, lty = 1,
                    xlab = "threshold quantile", ylab = paste(x$criterion, "criterion"),
                    col = grDevices::hcl.colors(max(length(lam), 2L), "Zissou 1"),
                    ...)
  graphics::points(x$probs, x$best$criterion, pch = 1, cex = 2.2, lwd = 2)
  graphics::legend("topright", legend = paste("lambda =", format(lam, digits = 3)),
                   col = grDevices::hcl.colors(max(length(lam), 2L), "Zissou 1"),
                   lty = 1, pch = 16, bty = "n", cex = 0.8)
  invisible(x)
}
