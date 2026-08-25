test_that("greedy_peeling matches a literal R implementation", {
  set.seed(7)
  for (trial in 1:15) {
    p <- sample(15:60, 1)
    A <- matrix(rnorm(p * p), p)
    W <- (A + t(A)) / 2
    diag(W) <- 0
    thr <- stats::quantile(W[upper.tri(W)], runif(1, 0.4, 0.9))
    lam <- runif(1, 0, 1)

    Wt <- W
    Wt[Wt < thr] <- 0
    diag(Wt) <- 0

    got <- greedy_peeling(W, lambda = lam, threshold = thr)
    want <- ref_peel(Wt, lam)

    expect_identical(got$keep, want$keep)
    expect_identical(got$drop, want$drop)
  }
})

test_that("greedy_peeling partitions the nodes exactly once", {
  set.seed(11)
  W <- cor(matrix(rnorm(60 * 40), 60))
  res <- greedy_peeling(W, lambda = 0.5, threshold = 0.1)
  expect_setequal(c(res$keep, res$drop), seq_len(40))
  expect_length(intersect(res$keep, res$drop), 0)
})

test_that("icons_detect returns a valid partition of every variable", {
  data(sim)
  W <- cor(sim)
  for (thr in c(0.3, 0.5, 0.7)) {
    part <- icons_detect(W, threshold = thr)
    expect_s3_class(part, "icons_partition")
    # order is a permutation
    expect_setequal(part$order, seq_len(ncol(W)))
    expect_length(part$order, ncol(W))
    # sizes and singletons account for everything
    expect_equal(sum(part$sizes) + part$n_singletons, ncol(W))
    # membership agrees with order and sizes
    expect_equal(sum(part$membership > 0), sum(part$sizes))
    expect_equal(unname(tabulate(part$membership, nbins = length(part$sizes))),
                 part$sizes)
    # min_size is respected
    expect_true(all(part$sizes >= 2))
  }
})

test_that("icons_detect is deterministic", {
  data(sim)
  W <- cor(sim)
  expect_identical(icons_detect(W, threshold = 0.5)$order,
                   icons_detect(W, threshold = 0.5)$order)
})

test_that("probs and threshold agree", {
  data(sim)
  W <- cor(sim)
  thr <- unname(stats::quantile(W[upper.tri(W)], 0.95))
  expect_equal(icons_detect(W, probs = 0.95)$sizes,
               icons_detect(W, threshold = thr)$sizes)
})

test_that("min_size and max_k are honoured", {
  data(sim)
  W <- cor(sim)
  expect_true(all(icons_detect(W, threshold = 0.5, min_size = 8)$sizes >= 8))
  expect_lte(length(icons_detect(W, threshold = 0.5, max_k = 3)$sizes), 3)
})

test_that("larger lambda gives smaller communities", {
  data(sim)
  W <- cor(sim)
  lo <- icons_detect(W, threshold = 0.4, lambda = 0.4)$sizes[1]
  hi <- icons_detect(W, threshold = 0.4, lambda = 0.9)$sizes[1]
  expect_gte(lo, hi)
})

test_that("input validation catches malformed matrices", {
  expect_error(icons_detect(matrix(1:6, 2)), class = "icons_dim_error")
  W <- matrix(1:9, 3)
  expect_error(icons_detect(W), class = "icons_sym_error")
  W2 <- diag(3); W2[1, 2] <- W2[2, 1] <- NA
  expect_error(icons_detect(W2), class = "icons_na_error")
  expect_error(icons_detect(diag(3), lambda = 2), class = "icons_value_error")
  expect_error(icons_detect(diag(3), min_size = 1), class = "icons_value_error")
  expect_error(icons_detect("a"), class = "icons_type_error")
})

test_that("reorder_matrix and block_index are consistent", {
  data(sim)
  W <- cor(sim)
  part <- icons_detect(W, threshold = 0.6)
  R <- reorder_matrix(W, part)
  expect_equal(dim(R), dim(W))
  expect_equal(R, W[part$order, part$order])

  bi <- block_index(part, 1)
  expect_length(bi$index, part$sizes[1])
  expect_true(all(part$membership[bi$index] == 1L))

  full <- reorder_matrix(W, part, communities_only = TRUE)
  expect_equal(ncol(full), sum(part$sizes))
})

test_that("as_membership round-trips", {
  data(sim)
  part <- icons_detect(cor(sim), threshold = 0.6)
  m <- as_membership(part)
  expect_length(m, part$p)
  expect_equal(unname(m), part$membership)
  expect_true(all(is.na(as_membership(part, singleton_na = TRUE)[m == 0])))
})
