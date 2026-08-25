#' Choose the number of factors
#'
#' Traces the fit criterion as communities are admitted to the model one at a
#' time, largest first, and finds the elbow.  With \eqn{k = 0} no variable
#' loads on a factor and the criterion equals \eqn{\|\mathrm{offdiag}(S)\|_F};
#' each additional community explains more covariance, with diminishing
#' returns once the real communities are used up.
#'
#' @section Cost:
#' The whole path costs one \eqn{O(n^2p)} pass plus \eqn{O(K^3)}, not \eqn{K}
#' separate fits.  The block summaries that the estimators depend on are
#' computed once for the full partition, and the sub-model for \eqn{k}
#' communities reuses the leading sub-blocks unchanged --- \eqn{\hat a_{kk}}
#' and \eqn{\hat b_{kk'}} depend only on their own blocks, so admitting further
#' communities never alters them.
#'
#' @param data A numeric matrix or data frame, observations in rows.
#' @param partition An `"icons_partition"` from [icons_detect()], or a
#'   membership vector.
#' @param k_max Largest number of factors to consider.  Defaults to every
#'   community in `partition`.
#' @param criterion `"frobenius"` (default) uses the fast incremental path;
#'   `"legacy"` refits the 0.1.x objective at every `k` and is much slower.
#'
#' @return An object of class `"icons_nfactors"`:
#'   \describe{
#'     \item{`path`}{Data frame with `k`, `criterion`, `relative`, `p_modelled`
#'       and `n_singletons`.}
#'     \item{`suggested`}{Elbow of the curve: the `k` furthest from the chord
#'       joining its endpoints.}
#'   }
#'
#' @seealso [icons_detect()], [scfa()], [scfa_criterion()]
#'
#' @examples
#' data(sim)
#' part <- icons_detect(cor(sim), threshold = 0.6)
#' nf <- n_factors(sim, part)
#' nf
#' plot(nf)
#'
#' @export
n_factors <- function(data,
                      partition,
                      k_max = NULL,
                      criterion = c("frobenius", "legacy")) {
  cl <- match.call()
  criterion <- match.arg(criterion)
  X <- as_data_matrix(data, arg = "data")
  p <- ncol(X)
  memb <- resolve_membership(partition, p)
  K <- max(0L, max(memb))
  if (K < 1L) {
    abort_icons("`partition` contains no community.",
                class = "icons_value_error")
  }
  k_max <- if (is.null(k_max)) K else min(K, check_count(k_max, "k_max", 1L))

  if (criterion == "legacy") {
    path <- do.call(rbind, lapply(0:k_max, function(k) {
      if (k == 0L) {
        st <- scfa_stats(X, integer(p), 0L)
        n2S <- norm2_S(X, st$center, st$den)
        val <- sqrt(max(n2S - sum(st$svar^2), 0))
        return(data.frame(k = 0L, criterion = val, relative = val / sqrt(n2S),
                          p_modelled = 0L, n_singletons = p))
      }
      m <- memb; m[m > k] <- 0L
      fit <- suppressWarnings(scfa(X, m, criterion = "legacy"))
      data.frame(k = k, criterion = fit$criterion$value,
                 relative = fit$criterion$relative,
                 p_modelled = fit$p_modelled, n_singletons = fit$n_singletons)
    }))
  } else {
    path <- nfactors_fast(X, memb, K, k_max)
  }
  rownames(path) <- NULL

  structure(
    list(path = path, suggested = elbow_point(path$k, path$criterion),
         criterion = criterion, K = K, call = cl),
    class = "icons_nfactors"
  )
}

# Incremental criterion path: one pass over the data, then O(K^3).
nfactors_fast <- function(X, memb, K, k_max) {
  st <- scfa_stats(X, memb, K)
  par <- scfa_params(st)
  n2S <- norm2_S(X, st$center, st$den)

  ss2_all <- sum(st$svar^2)
  # Per-community sum of s_jj^2, so the singleton contribution is a difference.
  ss2k <- vapply(seq_len(K), function(k) sum(st$svar[memb == k]^2), numeric(1L))

  vals <- vapply(0:k_max, function(k) {
    if (k == 0L) return(sqrt(max(n2S - ss2_all, 0)))
    ix <- seq_len(k)
    sub <- list(sizes = st$sizes[ix], trk = st$trk[ix],
                Sblk = st$Sblk[ix, ix, drop = FALSE],
                svar = st$svar, membership = ifelse(memb %in% ix, memb, 0L))
    tt <- uniform_block_terms(sub, list(a = par$a[ix],
                                        B = par$B[ix, ix, drop = FALSE]))
    sqrt(max(n2S - 2 * tt$cross + tt$norm2, 0))
  }, numeric(1L))

  pm <- c(0L, cumsum(st$sizes)[seq_len(k_max)])
  data.frame(k = 0:k_max, criterion = vals, relative = vals / sqrt(n2S),
             p_modelled = pm, n_singletons = st$p - pm)
}

# Kneedle-style elbow: the point furthest from the chord joining the endpoints.
elbow_point <- function(k, y) {
  if (length(k) < 3L || anyNA(y)) return(NA_integer_)
  xs <- (k - min(k)) / diff(range(k))
  ys <- (y - min(y)) / diff(range(y))
  if (!is.finite(sum(ys))) return(NA_integer_)
  # Distance from the line through the first and last points.
  x1 <- xs[1L]; y1 <- ys[1L]
  x2 <- xs[length(xs)]; y2 <- ys[length(ys)]
  d <- abs((y2 - y1) * xs - (x2 - x1) * ys + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  as.integer(k[which.max(d)])
}

#' @export
print.icons_nfactors <- function(x, digits = 4L, ...) {
  cat("<ICONS factor-count selection>\n")
  cat("  communities available :", x$K, "\n")
  cat("  criterion             :", x$criterion, "\n")
  cat("  suggested k (elbow)   :", x$suggested, "\n\n")
  print(x$path, row.names = FALSE, digits = digits)
  invisible(x)
}

#' Plot the factor-count criterion path
#'
#' @param x An `"icons_nfactors"` object.
#' @param relative Logical; plot the criterion relative to \eqn{\|S\|_F}.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' data(sim)
#' nf <- n_factors(sim, icons_detect(cor(sim), threshold = 0.6))
#' plot(nf)
#'
#' @export
plot.icons_nfactors <- function(x, relative = FALSE, ...) {
  y <- if (relative) x$path$relative else x$path$criterion
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(x$path$k, y, type = "b", pch = 16,
                 xlab = "number of factors (k)",
                 ylab = if (relative) "relative criterion" else "criterion",
                 ...)
  if (!is.na(x$suggested)) {
    graphics::abline(v = x$suggested, lty = 2, col = "grey40")
    graphics::mtext(paste("elbow at k =", x$suggested), side = 3,
                    adj = 1, cex = 0.8, col = "grey40")
  }
  invisible(x)
}
