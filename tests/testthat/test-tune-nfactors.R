test_that("icons_tune searches the full grid and selects the minimum", {
  data(sim)
  W <- cor(sim)
  probs <- c(0.90, 0.94, 0.97)
  lambda <- c(0.5, 0.7)
  tuned <- icons_tune(W, sim, probs = probs, lambda = lambda)

  expect_s3_class(tuned, "icons_tune")
  expect_equal(nrow(tuned$grid), length(probs) * length(lambda))
  expect_setequal(tuned$grid$probs, probs)
  expect_setequal(tuned$grid$lambda, lambda)
  expect_equal(tuned$best$criterion, min(tuned$grid$criterion, na.rm = TRUE))
  expect_true(tuned$threshold %in% tuned$grid$threshold)

  # The selected settings reproduce the reported structure.
  part <- icons_detect(W, threshold = tuned$threshold, lambda = tuned$lambda)
  expect_equal(length(part$sizes), tuned$best$k)
})

test_that("a grid cell that yields no community is recorded, not fatal", {
  set.seed(31)
  # Pure noise: at a high threshold nothing survives.
  W <- cor(matrix(rnorm(40 * 60), 40))
  X <- matrix(rnorm(40 * 60), 40)
  expect_error(
    tuned <- icons_tune(W, X, probs = c(0.5, 0.999), lambda = 0.5),
    NA
  )
  expect_equal(nrow(tuned$grid), 2L)
})

test_that("icons_tune validates its inputs", {
  data(sim)
  W <- cor(sim)
  expect_error(icons_tune(W, sim[, 1:10]), class = "icons_dim_error")
  expect_error(icons_tune(W, sim, probs = 1.5), class = "icons_value_error")
  expect_error(icons_tune(W, sim, lambda = -1), class = "icons_value_error")
})

test_that("n_factors traces a decreasing criterion and finds an elbow", {
  data(sim)
  part <- icons_detect(cor(sim), threshold = 0.6)
  nf <- n_factors(sim, part)

  expect_s3_class(nf, "icons_nfactors")
  expect_equal(nf$path$k, 0:length(part$sizes))
  # Adding a community can only explain more covariance.
  expect_true(all(diff(nf$path$criterion) <= 1e-8))
  # k = 0 is the off-diagonal norm of S.
  S <- stats::cov(sim)
  diag(S) <- 0
  expect_equal(nf$path$criterion[1], sqrt(sum(S^2)))
  expect_true(nf$suggested >= 1 && nf$suggested <= length(part$sizes))
  expect_equal(nf$path$p_modelled, c(0, cumsum(part$sizes)))
})

test_that("the fast n_factors path equals refitting at each k", {
  data(sim)
  part <- icons_detect(cor(sim), threshold = 0.6)
  nf <- n_factors(sim, part)

  slow <- vapply(seq_along(part$sizes), function(k) {
    m <- part$membership
    m[m > k] <- 0L
    suppressWarnings(scfa(sim, m)$criterion$value)
  }, numeric(1))

  expect_equal(nf$path$criterion[-1], slow)
})

test_that("n_factors respects k_max", {
  data(sim)
  part <- icons_detect(cor(sim), threshold = 0.6)
  expect_equal(max(n_factors(sim, part, k_max = 3)$path$k), 3L)
})

test_that("the legacy criterion reproduces the ICONS 0.1.x objective", {
  # ICONS 0.1.9 computed  ||offdiag(cov(X - F L'))||_F +
  #                       ||offdiag(cov(X) - L diag(cov(F)) L')||_F
  # with F the precision-weighted community means.  Recreate it literally.
  set.seed(37)
  memb <- rep(1:3, c(8, 6, 9))
  p <- length(memb)
  X <- matrix(rnorm(50 * p), 50)
  for (k in 1:3) X[, memb == k] <- X[, memb == k] + rnorm(50)

  fit <- scfa(X, memb, criterion = "legacy")

  Xc <- sweep(X, 2, colMeans(X))
  L <- factor_loadings(fit)
  Fols <- Xc %*% L %*% diag(1 / fit$sizes)
  psi <- colMeans((Xc - Fols[, memb])^2) * nrow(X) / (nrow(X) - 1)
  Fw <- sapply(1:3, function(k) {
    w <- 1 / psi[memb == k]
    as.numeric(Xc[, memb == k] %*% w) / sum(w)
  })

  U <- Xc - Fw[, memb]
  Cu <- stats::cov(U); diag(Cu) <- 0
  term1 <- sqrt(sum(Cu^2))

  Sf0 <- diag(diag(stats::cov(Fw)))
  D <- stats::cov(X) - L %*% Sf0 %*% t(L)
  diag(D) <- 0
  term2 <- sqrt(sum(D^2))

  expect_equal(fit$criterion$value, term1 + term2)
})

test_that("deprecated functions still work and warn", {
  data(sim)
  W <- cor(sim)
  expect_warning(old <- dense(W, 0.6, 0.5), "deprecated")
  expect_equal(sum(old$CID), ncol(W))
  expect_setequal(old$Clist, seq_len(ncol(W)))
  expect_equal(dim(old$W_dense), dim(W))

  expect_warning(m <- get_membership(c(2, 3), c(1, 4, 2, 5, 7)), "deprecated")
  expect_equal(m, c(1L, 2L, 0L, 1L, 2L, 0L, 2L))

  expect_warning(v <- get_vectorform(W), "deprecated")
  expect_equal(v, W[upper.tri(W)])

  expect_warning(gi <- get_index(1:2, c(3, 4, 2), c(3, 1, 2, 6, 4, 5, 7, 9, 8)),
                 "deprecated")
  expect_equal(gi$indices, c(3, 1, 2, 6, 4, 5, 7))
})
