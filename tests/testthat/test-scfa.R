test_that("estimators match equation (4) computed from the full covariance", {
  set.seed(11)
  n <- 90
  memb <- sample(c(rep(1:4, c(9, 6, 11, 4)), rep(0L, 7)))
  p <- length(memb)
  X <- matrix(rnorm(n * p), n)
  for (k in 1:4) X[, memb == k] <- X[, memb == k] + rnorm(n) * sqrt(1 + k)

  fit <- scfa(X, memb)
  want <- ref_scfa_params(X, memb)

  expect_equal(unname(fit$a), want$a)
  expect_equal(unname(fit$Sigma_f), want$B)
})

test_that("factor scores are the community means and equal the GLS solution", {
  set.seed(3)
  sim <- sim_scfa(80, c(6, 9, 5), c(0.2, 0.4, 0.3),
                  diag(3) + 0.3)
  fit <- scfa(sim$X, sim$membership)

  L <- factor_loadings(fit)
  Xc <- sweep(sim$X, 2, colMeans(sim$X))

  ols <- Xc %*% L %*% diag(1 / fit$sizes)
  expect_equal(unname(fit$scores), ols)

  Su <- diag(fit$a[sim$membership])
  gls <- t(solve(t(L) %*% solve(Su) %*% L, t(L) %*% solve(Su) %*% t(Xc)))
  expect_equal(unname(fit$scores), unname(gls))
})

test_that("cov(F_hat) equals Sigma_f + diag(a / p_k) exactly (Theorem 3.2)", {
  set.seed(5)
  sim <- sim_scfa(120, c(8, 12, 6), c(0.1, 0.2, 0.5), diag(3) * 2 + 0.5)
  fit <- scfa(sim$X, sim$membership)
  expect_equal(unname(stats::cov(fit$scores)),
               unname(fit$Sigma_f + diag(fit$a / fit$sizes)))
  # Which is why cov(F_hat), as used by ICONS 0.1.x, overstates the diagonal.
  expect_true(all(diag(stats::cov(fit$scores)) > diag(fit$Sigma_f)))
})

test_that("estimators are unbiased and Wald intervals cover at the nominal rate", {
  skip_on_cran()
  set.seed(2024)
  sizes <- c(10, 12, 8)
  a <- c(0.1, 0.2, 0.5)
  B <- matrix(c(2.02, 0.73, 1.15,
                0.73, 3.13, 1.63,
                1.15, 1.63, 3.69), 3, 3)
  n <- 120
  reps <- 300

  ahat <- matrix(NA_real_, reps, 3)
  cover_a <- matrix(NA, reps, 3)
  cover_b <- matrix(NA, reps, 3)
  for (r in seq_len(reps)) {
    d <- sim_scfa(n, sizes, a, B)
    f <- scfa(d$X, d$membership, criterion = "none")
    ci <- confint(f)
    ahat[r, ] <- f$a
    cover_a[r, ] <- a >= ci[1:3, 3] & a <= ci[1:3, 4]
    bd <- c("b[F1,F1]", "b[F2,F2]", "b[F3,F3]")
    cover_b[r, ] <- diag(B) >= ci[bd, 3] & diag(B) <= ci[bd, 4]
  }

  # Unbiasedness: the Monte Carlo mean is within 4 standard errors of truth.
  mc_se <- apply(ahat, 2, sd) / sqrt(reps)
  expect_true(all(abs(colMeans(ahat) - a) < 4 * mc_se))

  # Nominal 95% coverage, allowing Monte Carlo error.
  expect_true(all(abs(colMeans(cover_a) - 0.95) < 0.05))
  expect_true(all(abs(colMeans(cover_b) - 0.95) < 0.06))
})

test_that("the Frobenius criterion matches an explicit p by p computation", {
  set.seed(13)
  memb <- sample(c(rep(1:3, c(7, 9, 5)), rep(0L, 6)))
  X <- matrix(rnorm(70 * length(memb)), 70)
  fit <- scfa(X, memb)
  want <- ref_frobenius(X, memb, unname(fit$a), unname(fit$Sigma_f))
  expect_equal(fit$criterion$value, want)
  expect_equal(fit$criterion$relative, want / fit$criterion$norm_S)
})

test_that("criterion = 'none' skips the computation", {
  data(sim)
  fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6), criterion = "none")
  expect_true(is.na(scfa_criterion(fit)))
})

test_that("sigma_u = 'diagonal' gives precision-weighted scores", {
  set.seed(17)
  memb <- rep(1:3, c(8, 6, 10))
  X <- matrix(rnorm(60 * length(memb)), 60)
  X[, memb == 1] <- X[, memb == 1] * rep(c(1, 5), length.out = 8)  # heteroscedastic
  fit <- scfa(X, memb, sigma_u = "diagonal")

  psi <- sigma_u(fit)
  expect_length(psi, ncol(X))
  Xc <- sweep(X, 2, colMeans(X))
  for (k in 1:3) {
    w <- 1 / psi[memb == k]
    expect_equal(unname(fit$scores[, k]),
                 as.numeric(Xc[, memb == k] %*% w) / sum(w))
  }
  # Under the MLE the weights are constant within a community, so the two
  # score estimators coincide there.
  fit2 <- scfa(X, memb)
  expect_equal(unname(fit2$scores),
               unname(sweep(Xc %*% factor_loadings(fit2), 2, fit2$sizes, "/")))
})

test_that("partition objects and membership vectors agree", {
  data(sim)
  part <- icons_detect(cor(sim), threshold = 0.6)
  expect_equal(scfa(sim, part)$a, scfa(sim, part$membership)$a)
})

test_that("membership labels may be arbitrary and are relabelled", {
  set.seed(19)
  memb <- rep(c(50L, 7L, 0L), c(6, 8, 4))
  X <- matrix(rnorm(40 * 18), 40)
  fit <- scfa(X, memb)
  expect_equal(fit$K, 2L)
  expect_equal(fit$sizes, c(8L, 6L))   # community 7 relabelled to 1
  expect_equal(fit$n_singletons, 4L)
})

test_that("invalid inputs are rejected", {
  X <- matrix(rnorm(100), 20)
  # wrong length
  expect_error(scfa(X, rep(1L, 3)), class = "icons_type_error")
  # no community at all
  expect_error(scfa(X, rep(0L, 5)), class = "icons_value_error")
  # a community of one variable: a_kk would divide by zero
  expect_error(scfa(X, c(1L, 1L, 2L, 0L, 0L)), class = "icons_value_error")
  # negative labels
  expect_error(scfa(X, c(1L, 1L, -1L, 0L, 0L)), class = "icons_value_error")
  Xna <- X; Xna[1, 1] <- NA
  expect_error(scfa(Xna, rep(1:2, c(3, 2))), class = "icons_na_error")
})

test_that("q >= n warns", {
  set.seed(23)
  memb <- rep(1:6, each = 3)
  X <- matrix(rnorm(15 * 18), 15)
  expect_warning(scfa(X, memb), "q < n")
})
