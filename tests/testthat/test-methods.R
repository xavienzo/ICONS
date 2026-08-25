make_fit <- function() {
  data(sim, package = "ICONS", envir = environment())
  scfa(sim, icons_detect(cor(sim), threshold = 0.6))
}

test_that("coef, vcov and confint are mutually consistent", {
  fit <- make_fit()
  est <- coef(fit)
  V <- vcov(fit)
  ci <- confint(fit, level = 0.95)

  expect_length(est, fit$K + fit$K * (fit$K + 1) / 2)
  expect_equal(rownames(V), names(est))
  expect_equal(ci[, "estimate"], est)
  expect_equal(unname(ci[, "se"]), unname(sqrt(diag(V))))
  z <- qnorm(0.975)
  expect_equal(unname(ci[, 3]), unname(est - z * sqrt(diag(V))))
  expect_equal(unname(ci[, 4]), unname(est + z * sqrt(diag(V))))

  expect_equal(coef(fit, "a"), fit$a)
  expect_equal(coef(fit, "Sigma_f"), fit$Sigma_f)
  expect_equal(nrow(confint(fit, parm = 1:3)), 3L)
})

test_that("standard errors match Theorem 3 term by term", {
  set.seed(29)
  sizes <- c(9, 7, 11)
  a <- c(0.3, 0.6, 0.2)
  B <- matrix(c(2, 0.5, 0.4, 0.5, 3, 0.7, 0.4, 0.7, 1.5), 3, 3)
  d <- sim_scfa(100, sizes, a, B)
  fit <- scfa(d$X, d$membership, criterion = "none")

  n <- 100
  pk <- fit$sizes
  ah <- unname(fit$a)
  Bh <- unname(fit$Sigma_f)

  expect_equal(unname(fit$se$a), sqrt(2 * ah^2 / ((n - 1) * (pk - 1))))

  g <- ah + pk * diag(Bh)
  expect_equal(unname(diag(fit$se$B)),
               sqrt(2 / ((n - 1) * pk * (pk - 1)) *
                      (g^2 - (2 * ah + pk * diag(Bh)) * diag(Bh))))

  k <- 1; l <- 3
  want <- sqrt((pk[k] * pk[l] * (Bh[k, l]^2 + Bh[l, k]^2) +
                  2 * g[k] * g[l]) / (2 * (n - 1) * pk[k] * pk[l]))
  expect_equal(unname(fit$se$B[k, l]), want)
})

test_that("cov_scores is Sigma_f + diag(a / p_k)", {
  fit <- make_fit()
  expect_equal(fit$cov_scores, fit$Sigma_f + diag(fit$a / fit$sizes))
})

test_that("factor_loadings builds a valid indicator matrix", {
  fit <- make_fit()
  L <- factor_loadings(fit)
  expect_equal(dim(L), c(fit$p, fit$K))
  expect_true(all(L %in% c(0, 1)))
  expect_equal(unname(colSums(L)), fit$sizes)
  expect_equal(sum(L), sum(fit$sizes))
  expect_true(all(rowSums(L) <= 1))          # non-overlapping loadings

  sp <- factor_loadings(fit, sparse = TRUE)
  expect_equal(nrow(sp), sum(fit$sizes))
})

test_that("predict reproduces the training scores and handles new data", {
  fit <- make_fit()
  data(sim, package = "ICONS", envir = environment())
  expect_equal(predict(fit), fit$scores)
  expect_equal(predict(fit, sim), fit$scores)
  expect_equal(dim(predict(fit, sim[1:5, ])), c(5L, fit$K))
  expect_error(predict(fit, sim[, 1:10]), class = "icons_dim_error")
})

test_that("residuals are centred and orthogonal to the fitted scores", {
  fit <- make_fit()
  data(sim, package = "ICONS", envir = environment())
  U <- residuals(fit, sim)
  expect_equal(dim(U), c(fit$n, fit$p))
  expect_true(max(abs(colMeans(U))) < 1e-10)
  # Within a community the residuals sum to zero across variables, because the
  # score is that community's mean.
  k1 <- which(fit$membership == 1L)
  expect_true(max(abs(rowSums(U[, k1]))) < 1e-10)
  expect_error(residuals(fit), class = "icons_value_error")
})

test_that("fitted reconstructs the model-implied covariance and guards on size", {
  fit <- make_fit()
  S <- fitted(fit)
  expect_equal(dim(S), c(fit$p, fit$p))
  expect_true(isSymmetric(unname(S)))
  i <- which(fit$membership == 1L)[1:2]
  expect_equal(unname(S[i[1], i[2]]), unname(fit$Sigma_f[1, 1]))
  expect_equal(unname(S[i[1], i[1]]), unname(fit$Sigma_f[1, 1] + fit$a[1]))
  expect_error(fitted(fit, max_p = 10), class = "icons_size_error")
})

test_that("sigma_u returns per-variable error variances", {
  fit <- make_fit()
  s <- sigma_u(fit)
  expect_length(s, fit$p)
  i <- which(fit$membership == 2L)
  expect_true(all(s[i] == fit$a[2]))
  expect_true(all(is.na(s[fit$membership == 0L])))
})

test_that("nobs and print/summary work", {
  fit <- make_fit()
  expect_equal(nobs(fit), fit$n)
  expect_output(print(fit), "Semi-confirmatory factor analysis")
  expect_output(print(summary(fit)), "Wald intervals")
  part <- icons_detect(cor(matrix(rnorm(3000), 30)), threshold = 0.3)
  expect_output(print(part), "ICONS partition")
  expect_output(print(summary(part)), "community")
})

test_that("half_vec returns the strict upper triangle", {
  W <- cor(matrix(rnorm(200), 20))
  v <- half_vec(W)
  expect_length(v, ncol(W) * (ncol(W) - 1) / 2)
  expect_equal(v, W[upper.tri(W)])
  expect_error(half_vec(matrix(1:6, 2)), class = "icons_dim_error")
})

test_that("plotting functions run and restore graphics state", {
  fit <- make_fit()
  data(sim, package = "ICONS", envir = environment())
  part <- icons_detect(cor(sim), threshold = 0.6)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  before <- par("mfrow")
  expect_silent(plot_matrix(cor(sim)))
  expect_equal(par("mfrow"), before)

  expect_silent(plot_matrix(reorder_matrix(cor(sim), part), partition = part))
  expect_silent(plot(part, cor(sim)))
  expect_silent(plot(n_factors(sim, part)))
  expect_error(plot_matrix(matrix("a", 2, 2)), class = "icons_type_error")
})

test_that("plot_matrix downsamples large matrices", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  W <- cor(matrix(rnorm(60 * 300), 60))
  expect_silent(plot_matrix(W, max_cells = 50))
})

test_that("axis labels are round numbers in original units, not rescaled ones", {
  # Regression: ticks used to be pretty() on the downsampled grid with the
  # labels multiplied back, which turned 4199 into 1050 / 2100 / 3149 / 4199.
  for (n in c(200L, 4199L, 1000L)) {
    at <- ICONS:::index_ticks(n)
    expect_true(all(at >= 1 & at <= n))
    expect_true(all(at %% 5 == 0))            # ends in 5 or 0
    expect_equal(length(unique(diff(at))), 1L)  # equally spaced
  }
  expect_equal(ICONS:::index_ticks(4199L), c(1000, 2000, 3000, 4000))
})

test_that("downsampling uses whole-number blocks and preserves the mean", {
  z <- matrix(1:400 + 0.5, 20, 20)
  d <- ICONS:::downsample(z, 4L)
  expect_equal(dim(d), c(5L, 5L))
  expect_equal(d[1, 1], mean(z[1:4, 1:4]))
  expect_equal(mean(d), mean(z))

  # A side that is not a multiple of the block size keeps a partial last cell.
  z2 <- matrix(1, 7, 7)
  d2 <- ICONS:::downsample(z2, 3L)
  expect_equal(dim(d2), c(3L, 3L))
  expect_true(all(d2 == 1))
})
