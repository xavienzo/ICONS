## Semi-confirmatory factor analysis ------------------------------------------
##
## Everything below rests on one observation.  The closed-form estimators of
## Yang, Ma, Bi and Chen (2024) depend on the sample covariance matrix S only
## through, for each pair of communities,
##
##     sum(S_kk')   and   tr(S_kk),
##
## and both follow from the within-community sums T = X_c L (n by K) and the
## per-variable sums of squares.  So S itself is never needed: the estimation
## cost is O(np + nK^2) rather than O(np^2), and no p by p object is ever
## allocated.  Sub-blocks of the same T and the same sums of squares serve
## every candidate partition, which is what makes icons_tune() and n_factors()
## affordable.

# Single pass over the data giving everything the estimators need.
# `membership` is in original variable order, 0 for singletons.
scfa_stats <- function(X, membership, K, center = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  mu <- if (center) .colMeans(X, n, p) else numeric(p)
  bs <- .block_stats_cpp(X, mu, as.integer(membership), as.integer(K))
  den <- if (center) n - 1 else n

  svar <- bs$ss / den                       # s_jj for every variable
  Tk <- bs$T                                # n by K, centred community sums
  Sblk <- crossprod(Tk) / den               # Sblk[k, l] = sum(S_kl)
  sizes <- tabulate(membership, nbins = K)
  # tr(S_kk), accumulated by community in one sweep.
  trk <- numeric(K)
  keep <- membership > 0L
  if (any(keep)) {
    agg <- rowsum(svar[keep], membership[keep], reorder = TRUE)
    trk[as.integer(rownames(agg))] <- as.numeric(agg)
  }

  list(n = n, p = p, K = K, sizes = sizes, svar = svar, Tk = Tk,
       Sblk = Sblk, trk = trk, center = mu, den = den,
       membership = membership)
}

# Closed-form MLEs (equation 4 of Yang, Ma, Bi and Chen, 2024) with tau = 1.
scfa_params <- function(st) {
  pk <- st$sizes
  d <- pk * (pk - 1)
  a <- (pk * st$trk - diag(st$Sblk)) / d
  B <- st$Sblk / outer(pk, pk)
  diag(B) <- (diag(st$Sblk) - st$trk) / d
  list(a = a, B = B)
}

# Exact variance estimators (Theorem 3 of Yang, Ma, Bi and Chen, 2024).
scfa_variances <- function(st, par) {
  n <- st$n
  pk <- st$sizes
  a <- par$a
  B <- par$B
  K <- st$K

  var_a <- 2 * a^2 / ((n - 1) * (pk - 1))

  var_B <- matrix(NA_real_, K, K, dimnames = dimnames(B))
  g <- a + pk * diag(B)                       # a_kk + p_k b_kk
  diag(var_B) <- 2 / ((n - 1) * pk * (pk - 1)) *
    (g^2 - (a + g) * diag(B))                 # 2a + p b = a + g
  if (K > 1L) {
    for (k in seq_len(K - 1L)) {
      for (l in seq(k + 1L, K)) {
        v <- (2 * pk[k] * pk[l] * B[k, l]^2 + 2 * g[k] * g[l]) /
          (2 * (n - 1) * pk[k] * pk[l])
        var_B[k, l] <- var_B[l, k] <- v
      }
    }
  }
  list(a = var_a, B = var_B)
}


#' Semi-confirmatory factor analysis
#'
#' Fits the SCFA model of Yang, Ma, Bi and Chen (2024).  The communities in
#' `partition` specify which loadings of \eqn{L} are non-zero, one common
#' factor per community, so the "confirmatory" structure is learned from the
#' data rather than assumed.  Every parameter then has a closed-form maximum
#' likelihood estimator, so no iterative optimiser is involved.
#'
#' @section Model:
#' With variables ordered by community,
#' \deqn{X = Lf + u, \qquad \Sigma = L\Sigma_f L^\top + \Sigma_u,}
#' \eqn{L = \mathrm{Bdiag}(1_{p_1}, \ldots, 1_{p_K})},
#' \eqn{\Sigma_f = (b_{kk'})} and
#' \eqn{\Sigma_u = \mathrm{Bdiag}(a_{11}I_{p_1}, \ldots, a_{KK}I_{p_K})}.
#' Writing \eqn{S_{kk'}} for the corresponding block of the sample covariance
#' matrix, the maximum likelihood estimators are
#' \deqn{\hat a_{kk} = \frac{p_k\,\mathrm{tr}(S_{kk}) - \mathrm{sum}(S_{kk})}{p_k(p_k-1)},
#'       \qquad
#'       \hat b_{kk} = \frac{\mathrm{sum}(S_{kk}) - \mathrm{tr}(S_{kk})}{p_k(p_k-1)},
#'       \qquad
#'       \hat b_{kk'} = \frac{\mathrm{sum}(S_{kk'})}{p_k p_{k'}} \ (k \neq k').}
#' These are uniformly minimum-variance unbiased whenever
#' \eqn{K + K(K+1)/2 < n}, and the factor scores
#' \eqn{\hat f_i = (L^\top L)^{-1} L^\top X_i} --- simply the community means of
#' the centred observations --- coincide with the GLS and feasible GLS
#' estimators.
#'
#' @section Singletons:
#' Variables that `partition` did not assign to any community are carried
#' through but not modelled: they load on no factor, and their variance enters
#' the fit criterion as wholly unexplained covariance.  Communities must have
#' at least two variables, since \eqn{\hat a_{kk}} divides by
#' \eqn{p_k(p_k-1)}.
#'
#' @param data A numeric matrix or data frame, observations in rows, variables
#'   in columns.  Variables are used in their original order; `partition`
#'   supplies the community assignment.
#' @param partition An `"icons_partition"` from [icons_detect()], or an integer
#'   membership vector of length `ncol(data)` using `0` or `NA` for
#'   unassigned variables.
#' @param sigma_u Error-variance model.  `"mle"` (default) is the closed-form
#'   UMVUE \eqn{\hat a_{kk}I_{p_k}}, constant within a community.  `"diagonal"`
#'   estimates an unconstrained variance per variable from the residuals; it is
#'   useful for checking the within-community homogeneity assumption, but it is
#'   not the UMVUE and it breaks the OLS/GLS equivalence for factor scores.
#' @param center Logical; centre the columns before fitting.  The model assumes
#'   a zero mean, so leave this `TRUE` unless the data are already centred.
#' @param criterion Fit criterion recorded with the model, see
#'   [scfa_criterion()].  `"frobenius"` (default) is
#'   \eqn{\|S - \hat\Sigma\|_F}; `"legacy"` reproduces the two-term objective
#'   used by ICONS 0.1.x.  `"none"` skips it, which avoids the \eqn{O(n^2p)}
#'   Gram matrix.
#'
#' @return An object of class `"scfa"`:
#'   \describe{
#'     \item{`a`}{Error variances \eqn{\hat a_{kk}}, one per community.}
#'     \item{`Sigma_f`}{\eqn{K \times K} estimated factor covariance matrix
#'       \eqn{(\hat b_{kk'})}.}
#'     \item{`scores`}{\eqn{n \times K} matrix of estimated factor scores.}
#'     \item{`se`}{Exact standard errors of `a` and `Sigma_f`.}
#'     \item{`cov_scores`}{Exact covariance matrix of a factor score vector,
#'       \eqn{\mathrm{diag}(a_{kk}/p_k) + \Sigma_f}.}
#'     \item{`criterion`}{Value of the recorded fit criterion, and its
#'       relative version \eqn{\|S - \hat\Sigma\|_F / \|S\|_F}.}
#'     \item{`sizes`, `membership`, `K`, `n`, `p`, `n_singletons`}{Structure of
#'       the fit.}
#'   }
#'   Loadings are not stored: \eqn{\hat L} is fully determined by the
#'   membership, and [factor_loadings()] builds it on demand.
#'
#' @references
#' Yang, Y., Ma, T., Bi, C., & Chen, S. (2024). Semi-confirmatory factor
#' analysis for high-dimensional data with interconnected community structures.
#' *arXiv*. \doi{10.48550/arXiv.2401.00624}
#'
#' Yang, Y., Chen, C., & Chen, S. (2024). Covariance matrix estimation for
#' high-throughput biomedical data with interconnected communities.
#' *The American Statistician*, 78(4), 401--411.
#' \doi{10.1080/00031305.2024.2329681}
#'
#' @seealso [icons_detect()], [n_factors()], [confint.scfa()], [factor_loadings()]
#'
#' @examples
#' data(sim)
#' part <- icons_detect(cor(sim), threshold = 0.6)
#' fit <- scfa(sim, part)
#' fit
#'
#' # Wald confidence intervals for the covariance parameters
#' head(confint(fit))
#'
#' # Factor scores
#' dim(fit$scores)
#'
#' @export
scfa <- function(data,
                 partition,
                 sigma_u = c("mle", "diagonal"),
                 center = TRUE,
                 criterion = c("frobenius", "legacy", "none")) {
  cl <- match.call()
  sigma_u <- match.arg(sigma_u)
  criterion <- match.arg(criterion)

  X <- as_data_matrix(data, arg = "data")
  n <- nrow(X)
  p <- ncol(X)

  memb <- resolve_membership(partition, p)
  K <- max(0L, max(memb))
  if (K < 1L) {
    abort_icons("`partition` contains no community; nothing to fit.",
                class = "icons_value_error")
  }
  sizes <- tabulate(memb, nbins = K)
  if (any(sizes < 2L)) {
    bad <- which(sizes < 2L)
    abort_icons("every community needs at least 2 variables; ",
                "community/communities ", paste(bad, collapse = ", "),
                " have fewer. Increase `min_size` in icons_detect().",
                class = "icons_value_error")
  }
  q <- K + K * (K + 1L) / 2L
  if (q >= n) {
    warning("K + K(K+1)/2 = ", q, " is not smaller than n = ", n,
            "; the UMVUE and Wald results of Yang et al. (2024) need q < n.",
            call. = FALSE)
  }

  st <- scfa_stats(X, memb, K, center = center)
  par <- scfa_params(st)

  if (any(par$a <= 0)) {
    warning("estimated error variance is non-positive for community/communities ",
            paste(which(par$a <= 0), collapse = ", "),
            "; Sigma is not positive definite and the fit is unreliable.",
            call. = FALSE)
  }

  # OLS scores: the community means.  Under sigma_u = "mle" these are also the
  # GLS and FGLS estimators (Theorem 2), so no reweighting is needed.
  scores <- sweep(st$Tk, 2L, st$sizes, "/")

  # An unconstrained diagonal Sigma_u breaks that equivalence, so the scores
  # must actually be reweighted: one feasible GLS step, giving precision-
  # weighted community means.  The legacy criterion needs the same quantity
  # because that is what ICONS 0.1.x reported.
  psi <- NULL
  need_weighted <- sigma_u == "diagonal" || criterion == "legacy"
  if (need_weighted) {
    psi <- residual_var(X, st, scores, memb)
    scores_w <- weighted_scores(X, st, psi, memb, K)
    if (sigma_u == "diagonal") scores <- scores_w
  }

  dimnames(scores) <- list(rownames(X), paste0("F", seq_len(K)))
  names(par$a) <- paste0("F", seq_len(K))
  dimnames(par$B) <- list(paste0("F", seq_len(K)), paste0("F", seq_len(K)))

  se <- scfa_variances(st, par)
  se <- list(a = sqrt(pmax(se$a, 0)), B = sqrt(pmax(se$B, 0)))
  names(se$a) <- names(par$a)
  dimnames(se$B) <- dimnames(par$B)

  if (sigma_u == "diagonal") {
    names(psi) <- colnames(X)
  } else {
    psi <- NULL
  }

  crit <- if (criterion == "none") {
    list(criterion = criterion, value = NA_real_, relative = NA_real_,
         norm_S = NA_real_)
  } else if (criterion == "legacy") {
    scfa_criterion_internal(X, st, par, scores_w, memb, criterion)
  } else {
    scfa_criterion_internal(X, st, par, scores, memb, criterion)
  }

  cov_scores <- par$B + diag(par$a / st$sizes, K)
  dimnames(cov_scores) <- dimnames(par$B)

  structure(
    list(
      a            = par$a,
      Sigma_f      = par$B,
      Sigma_u_diag = psi,
      scores       = scores,
      se           = se,
      cov_scores   = cov_scores,
      criterion    = crit,
      sizes        = st$sizes,
      membership   = memb,
      K            = K,
      n            = n,
      p            = p,
      p_modelled   = sum(st$sizes),
      n_singletons = sum(memb == 0L),
      sigma_u      = sigma_u,
      center       = st$center,
      var_names    = colnames(X),
      call         = cl
    ),
    class = "scfa"
  )
}


# Accept either an icons_partition or a bare membership vector.
resolve_membership <- function(partition, p) {
  if (inherits(partition, "icons_partition")) {
    if (partition$p != p) {
      abort_icons("`partition` describes ", partition$p,
                  " variables but `data` has ", p, " columns.",
                  class = "icons_dim_error")
    }
    return(partition$membership)
  }
  if (!is.numeric(partition) || length(partition) != p) {
    abort_icons("`partition` must be an <icons_partition> or an integer ",
                "vector of length ", p, ".", class = "icons_type_error")
  }
  m <- as.integer(partition)
  m[is.na(m)] <- 0L
  if (any(m < 0L)) {
    abort_icons("`partition` must not contain negative labels.",
                class = "icons_value_error")
  }
  # Relabel so communities are 1..K with no gaps.
  lab <- sort(unique(m[m > 0L]))
  out <- integer(p)
  out[m > 0L] <- match(m[m > 0L], lab)
  out
}

# Per-variable residual variance, var(X_j - f_hat_{phi(j)}).
residual_var <- function(X, st, scores, memb) {
  n <- st$n
  p <- st$p
  out <- numeric(p)
  for (j in seq_len(p)) {
    z <- X[, j] - st$center[j]
    k <- memb[j]
    r <- if (k > 0L) z - scores[, k] else z
    out[j] <- sum(r * r) / (n - 1)
  }
  # A constant variable gives psi = 0 and an infinite weight; floor it.
  pmax(out, .Machine$double.eps)
}

# Feasible GLS factor scores under an unconstrained diagonal Sigma_u:
# precision-weighted means of the community's variables.  Reduces to the plain
# community mean when the weights within a community are equal, which is
# exactly the sigma_u = "mle" case.
weighted_scores <- function(X, st, psi, memb, K) {
  out <- matrix(0, st$n, K)
  for (k in seq_len(K)) {
    idx <- which(memb == k)
    w <- 1 / psi[idx]
    Z <- sweep(X[, idx, drop = FALSE], 2L, st$center[idx],
               check.margin = FALSE)
    out[, k] <- as.numeric(Z %*% w) / sum(w)
  }
  out
}
