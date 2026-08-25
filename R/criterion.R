## Model-fit criteria ---------------------------------------------------------
##
## The criteria compare the sample covariance S with the model-implied
## covariance, both of which are p by p.  Neither is ever formed.  Two
## identities do the work:
##
##   * ||S||_F^2 = ||X_c' X_c||_F^2 / (n-1)^2 = ||X_c X_c'||_F^2 / (n-1)^2,
##     so the n by n Gram matrix suffices -- O(n^2 p) instead of O(n p^2), a
##     decisive difference when p >> n.
##
##   * the model-implied matrix is uniform-block, so its inner product with S
##     and its own squared norm collapse to sums over the K by K block
##     summaries already computed by scfa_stats().

# ||S||_F^2 from the Gram matrix of the centred data.
norm2_S <- function(X, center, den) {
  Xc <- sweep(X, 2L, center, check.margin = FALSE)
  G <- tcrossprod(Xc)
  sum(G * G) / den^2
}

# <S, M> and ||M||_F^2 for the model-implied covariance
#   M = L Sigma_f L' + Sigma_u   on the modelled variables,
#   M_jj = s_jj and M_jl = 0     on the singletons.
uniform_block_terms <- function(st, par) {
  pk <- st$sizes
  a <- par$a
  B <- par$B

  # <S, M>: within block k the diagonal of M is a_k + b_kk and the rest is
  # b_kk, hence a_k * tr(S_kk) + b_kk * sum(S_kk).
  cross <- sum(a * st$trk) + sum(B * st$Sblk)

  # ||M||_F^2: p_k entries equal to a_k + b_kk, p_k(p_k - 1) equal to b_kk,
  # and p_k p_l equal to b_kl off the diagonal.
  d <- diag(B)
  norm2 <- sum(pk * (a + d)^2 + pk * (pk - 1) * d^2) +
    (sum(outer(pk, pk) * B^2) - sum(pk^2 * d^2))

  # Singletons: M matches the diagonal exactly and is zero elsewhere.
  s_single <- st$svar[st$membership == 0L]
  cross <- cross + sum(s_single^2)
  norm2 <- norm2 + sum(s_single^2)

  list(cross = cross, norm2 = norm2)
}

scfa_criterion_internal <- function(X, st, par, scores, memb, criterion) {
  n2S <- norm2_S(X, st$center, st$den)

  if (criterion == "frobenius") {
    tt <- uniform_block_terms(st, par)
    val <- sqrt(max(n2S - 2 * tt$cross + tt$norm2, 0))
    return(list(criterion = criterion, value = val,
                relative = val / sqrt(n2S), norm_S = sqrt(n2S)))
  }

  ## Legacy objective of ICONS 0.1.x: the Frobenius norm of the off-diagonal
  ## residual covariance, plus the Frobenius norm of the off-diagonal part of
  ## S minus the within-community factor contribution.  Kept so that results
  ## published with 0.1.x can be reproduced; scfa_criterion() documents why
  ## "frobenius" is the better default.
  pk <- st$sizes
  Sf <- stats::cov(scores)              # what 0.1.x used in place of Sigma_f

  # Term 1: ||offdiag(cov(U))||_F with U = X_c - F_hat L'.
  U <- X
  U <- sweep(U, 2L, st$center, check.margin = FALSE)
  idx <- memb > 0L
  U[, idx] <- U[, idx] - scores[, memb[idx], drop = FALSE]
  Gu <- tcrossprod(U)
  n2U <- sum(Gu * Gu) / st$den^2
  du <- .colSums(U * U, st$n, st$p) / st$den
  term1 <- sqrt(max(n2U - sum(du^2), 0))

  # Term 2: ||offdiag(S - L diag(Sigma_f) L')||_F.
  d <- diag(Sf)
  n2D <- n2S - 2 * sum(d * diag(st$Sblk)) + sum(d^2 * pk^2)
  cj <- numeric(st$p)
  cj[idx] <- d[memb[idx]]
  term2 <- sqrt(max(n2D - sum((st$svar - cj)^2), 0))

  val <- term1 + term2
  list(criterion = criterion, value = val, relative = val / sqrt(n2S),
       norm_S = sqrt(n2S))
}


#' Fit criterion for an SCFA model
#'
#' How far the model-implied covariance sits from the sample covariance.
#'
#' `"frobenius"`, the default, is \eqn{\|S - \hat\Sigma\|_F} with
#' \eqn{\hat\Sigma = \hat L\hat\Sigma_f\hat L^\top + \hat\Sigma_u} on the
#' modelled variables and \eqn{\hat\Sigma_{jj} = s_{jj}} on the singletons.
#' It is the natural loss for the structured covariance model, it is on the
#' scale of \eqn{S}, and because unassigned variables contribute their whole
#' off-diagonal covariance it does not reward partitions that simply model
#' fewer variables.
#'
#' `"legacy"` reproduces the two-term objective minimised by
#' `param_tuning_sigmau()` in ICONS 0.1.x.  It adds two Frobenius norms that
#' measure overlapping quantities and uses `cov(F_hat)` in place of
#' \eqn{\hat\Sigma_f}, which Theorem 3 of Yang, Ma, Bi and Chen (2024) shows
#' overstates the factor variances by \eqn{a_{kk}/p_k}.  It is provided for
#' reproducing earlier results, not recommended for new work.
#'
#' Both are computed exactly, in \eqn{O(n^2p)} time and without forming any
#' \eqn{p \times p} matrix.
#'
#' @param object An `"scfa"` object.
#' @param relative Logical; return \eqn{\|S - \hat\Sigma\|_F / \|S\|_F} instead
#'   of the raw value, which is comparable across datasets.
#'
#' @return A single number.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' scfa_criterion(fit)
#' scfa_criterion(fit, relative = TRUE)
#'
#' @export
scfa_criterion <- function(object, relative = FALSE) {
  if (!inherits(object, "scfa")) {
    abort_icons("`object` must be an <scfa> fit, not ", class(object)[1L], ".",
                class = "icons_type_error")
  }
  if (relative) object$criterion$relative else object$criterion$value
}
