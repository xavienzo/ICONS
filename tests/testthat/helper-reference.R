# Deliberately naive reference implementations.  They mirror the published
# formulas as literally as possible -- forming S, looping over blocks -- so that
# the fast paths in the package are checked against the definitions rather than
# against themselves.

# Greedy peeling written straight from the description: strip the minimum-degree
# node, score the survivors, take the best cut.
ref_peel <- function(W, lambda) {
  N <- nrow(W)
  rec <- matrix(0, N, 2)
  C <- colSums(W)
  removed <- integer(0)
  for (ite in seq_len(N)) {
    C[removed] <- Inf
    j <- which.min(C)
    removed <- c(removed, j)
    C <- C - W[j, ]
    C[removed] <- 0
    rec[ite, ] <- c(j, (sum(C) / 2) / ((N - ite)^(2 * lambda)))
  }
  best <- which.max(rec[, 2])
  list(keep = as.integer(rec[N:(best + 1), 1]),
       drop = as.integer(rec[1:best, 1]))
}

# Equation (4) of Yang, Ma, Bi and Chen (2024), computed from the full p by p
# sample covariance matrix.
ref_scfa_params <- function(X, memb) {
  S <- stats::cov(X)
  K <- max(memb)
  a <- numeric(K)
  B <- matrix(0, K, K)
  for (k in seq_len(K)) {
    for (l in seq_len(K)) {
      pk <- sum(memb == k)
      pl <- sum(memb == l)
      Skl <- S[memb == k, memb == l, drop = FALSE]
      if (k == l) {
        a[k] <- (pk * sum(diag(Skl)) - sum(Skl)) / (pk * (pk - 1))
        B[k, k] <- (sum(Skl) - sum(diag(Skl))) / (pk * (pk - 1))
      } else {
        B[k, l] <- sum(Skl) / (pk * pl)
      }
    }
  }
  list(a = a, B = B, S = S)
}

# ||S - Sigma_hat||_F built explicitly as a p by p difference.
ref_frobenius <- function(X, memb, a, B) {
  S <- stats::cov(X)
  p <- ncol(X)
  idx <- memb > 0
  Sig <- matrix(0, p, p)
  Sig[idx, idx] <- B[memb[idx], memb[idx]]
  d <- numeric(p)
  d[idx] <- a[memb[idx]]
  d[!idx] <- diag(S)[!idx]
  diag(Sig) <- diag(Sig) + d
  sqrt(sum((S - Sig)^2))
}

# Draw from the SCFA generative model with known parameters.
sim_scfa <- function(n, sizes, a, B) {
  K <- length(sizes)
  p <- sum(sizes)
  g <- rep(seq_len(K), sizes)
  f <- matrix(stats::rnorm(n * K), n) %*% chol(B)
  u <- matrix(stats::rnorm(n * p), n) * rep(sqrt(a[g]), each = n)
  list(X = f %*% t(diag(K)[g, , drop = FALSE]) + u, membership = g)
}
