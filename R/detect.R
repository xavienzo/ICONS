#' Greedy peeling with \eqn{\ell_0} shrinkage
#'
#' Extracts a single dense subgraph.  Nodes are stripped one at a time, always
#' the one with the smallest remaining weighted degree, and the surviving set
#' is scored by the generalised density
#' \deqn{h_\lambda(S) = w(S) / |S|^{2\lambda},}
#' where \eqn{w(S)} is the total edge weight inside \eqn{S}.  The best-scoring
#' set over the whole peeling sequence is returned.  \eqn{\lambda = 1/2}
#' recovers the average-degree objective of Charikar (2000); larger
#' \eqn{\lambda} penalises size more heavily and yields smaller, denser
#' subgraphs.
#'
#' This is the inner step of [icons_detect()], which calls it repeatedly.  Use
#' `icons_detect()` unless you specifically want one subgraph.
#'
#' @param W A square symmetric numeric matrix of edge weights, typically a
#'   correlation or covariance matrix.  The diagonal is ignored.
#' @param lambda Shrinkage exponent \eqn{\lambda}, a single number in
#'   \eqn{[0, 1]}.  Defaults to `0.5`.
#' @param threshold Edge weights below `threshold` are set to zero before
#'   peeling.  Defaults to `-Inf`, which keeps every edge.
#'
#' @return A list with components
#'   \describe{
#'     \item{`keep`}{Integer node indices of the extracted subgraph, ordered
#'       from the densest core outwards.}
#'     \item{`drop`}{Integer node indices peeled away before the best cut.}
#'     \item{`score`}{The value of \eqn{h_\lambda} at the best cut.}
#'   }
#'
#' @references
#' Charikar, M. (2000). Greedy approximation algorithms for finding dense
#' components in a graph. *Approximation Algorithms for Combinatorial
#' Optimization*, 84--95.
#'
#' Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
#' Identifying covariate-related subnetworks for whole-brain connectome
#' analysis. *Biostatistics*, 25(2), 541--558.
#' \doi{10.1093/biostatistics/kxad007}
#'
#' Tsourakakis, C., Bonchi, F., Gionis, A., Gullo, F., & Tsiarli, M. (2013).
#' Denser than the densest subgraph: extracting optimal quasi-cliques with
#' quality guarantees. *KDD '13*, 104--112. \doi{10.1145/2487575.2487645}
#'
#' @seealso [icons_detect()]
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#' res <- greedy_peeling(W, lambda = 0.6, threshold = 0.5)
#' length(res$keep)
#'
#' @export
greedy_peeling <- function(W, lambda = 0.5, threshold = -Inf) {
  W <- as_weight_matrix(W, arg = "W")
  lambda <- check_scalar_number(lambda, "lambda", lower = 0, upper = 1)
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    abort_icons("`threshold` must be a single number.",
                class = "icons_type_error")
  }
  if (!is.finite(threshold)) threshold <- -.Machine$double.xmax
  .peel_cpp(W, lambda, as.numeric(threshold))
}


#' Detect interconnected communities in a covariance matrix
#'
#' Repeatedly applies [greedy_peeling()]: the densest subgraph is extracted and
#' removed, the algorithm runs again on what is left, and so on until the
#' densest remaining subgraph falls below `min_size`.  Variables never assigned
#' to a community form the *singleton set*, which is reported separately rather
#' than treated as another community.
#'
#' The resulting partition specifies the non-zero loadings of the
#' semi-confirmatory factor model fitted by [scfa()]; each community becomes one
#' common factor.
#'
#' `threshold` is the single most influential argument: it decides which edges
#' the peeling sees at all.  Supply it as a quantile of the observed edge
#' weights via `probs`, or tune it together with `lambda` using [icons_tune()].
#'
#' @param W A square symmetric numeric matrix, typically `cor(data)` or
#'   `cov(data)`.  The diagonal is ignored.
#' @param threshold Numeric edge threshold; weights below it are set to zero.
#'   Ignored when `probs` is supplied.  Defaults to `0.5`.
#' @param probs Optional single probability in \eqn{(0, 1)}.  When given, the
#'   threshold is the corresponding quantile of the strict upper triangle of
#'   `W`, e.g. `probs = 0.95` keeps the strongest 5% of edges.
#' @param lambda Shrinkage exponent \eqn{\lambda} in \eqn{[0, 1]}.  Defaults to
#'   `0.5`.
#' @param min_size Smallest community retained, at least `2`.  Anything smaller
#'   goes to the singleton set.  SCFA requires at least two variables per
#'   factor, so `2` is the lowest usable value.  Defaults to `2`.
#' @param max_k Optional cap on the number of communities.  `Inf` (the default)
#'   means no cap.
#'
#' @return An object of class `"icons_partition"`: a list with components
#'   \describe{
#'     \item{`sizes`}{Integer vector of community sizes, largest-first by
#'       extraction order.}
#'     \item{`order`}{Integer permutation of `seq_len(ncol(W))`: community
#'       members in block order, then the singleton set.}
#'     \item{`membership`}{Integer vector of length `ncol(W)` giving the
#'       community of each *original* variable, `0` for singletons.}
#'     \item{`n_singletons`}{Number of unassigned variables.}
#'     \item{`scores`}{Value of the peeling objective for each community.}
#'     \item{`threshold`, `lambda`, `min_size`}{The settings used.}
#'   }
#'
#' @references
#' Wu, Q., Huang, X., Culbreth, A. J., Waltz, J. A., Hong, L. E., & Chen, S.
#' (2022). Extracting brain disease-related connectome subgraphs by adaptive
#' dense subgraph discovery. *Biometrics*, 78(4), 1566--1578.
#' \doi{10.1111/biom.13537}
#'
#' Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
#' Identifying covariate-related subnetworks for whole-brain connectome
#' analysis. *Biostatistics*, 25(2), 541--558.
#' \doi{10.1093/biostatistics/kxad007}
#'
#' @seealso [icons_tune()] to choose `threshold` and `lambda`,
#'   [scfa()] to fit the factor model, [reorder_matrix()] and [plot_matrix()]
#'   to visualise the result.
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#'
#' part <- icons_detect(W, threshold = 0.6)
#' part
#'
#' # Threshold given as a quantile of the edge weights instead
#' icons_detect(W, probs = 0.95)
#'
#' @export
icons_detect <- function(W,
                         threshold = 0.5,
                         probs = NULL,
                         lambda = 0.5,
                         min_size = 2L,
                         max_k = Inf) {
  cl <- match.call()
  W <- as_weight_matrix(W, arg = "W")
  lambda <- check_scalar_number(lambda, "lambda", lower = 0, upper = 1)
  min_size <- check_count(min_size, "min_size", lower = 2L)

  if (!is.null(probs)) {
    probs <- check_scalar_number(probs, "probs", lower = 0, upper = 1)
    threshold <- unname(stats::quantile(W[upper.tri(W)], probs))
  } else {
    threshold <- check_scalar_number(threshold, "threshold")
  }

  if (is.infinite(max_k)) {
    max_k_i <- 0L
  } else {
    max_k_i <- check_count(max_k, "max_k", lower = 1L)
  }

  res <- .detect_cpp(W, threshold, lambda, min_size, max_k_i)

  p <- ncol(W)
  sizes <- as.integer(res$sizes)
  membership <- integer(p)
  if (length(sizes)) {
    membership[res$order[seq_len(sum(sizes))]] <- rep.int(seq_along(sizes), sizes)
  }

  structure(
    list(
      sizes        = sizes,
      order        = as.integer(res$order),
      membership   = membership,
      n_singletons = as.integer(res$n_singletons),
      scores       = as.numeric(res$scores),
      threshold    = threshold,
      lambda       = lambda,
      min_size     = min_size,
      p            = p,
      labels       = colnames(W),
      call         = cl
    ),
    class = "icons_partition"
  )
}


#' @export
print.icons_partition <- function(x, ...) {
  cat("<ICONS partition>\n")
  cat("  variables   :", x$p, "\n")
  cat("  communities :", length(x$sizes), "\n")
  cat("  singletons  :", x$n_singletons, "\n")
  cat("  threshold   :", format(x$threshold, digits = 4),
      "  lambda:", format(x$lambda, digits = 4), "\n")
  if (length(x$sizes)) {
    cat("  sizes       :", paste(utils::head(x$sizes, 15L), collapse = " "),
        if (length(x$sizes) > 15L) paste0("... (", length(x$sizes), " total)") else "",
        "\n")
  }
  invisible(x)
}

#' @export
summary.icons_partition <- function(object, ...) {
  k <- length(object$sizes)
  tab <- data.frame(
    community = seq_len(k),
    size      = object$sizes,
    score     = object$scores
  )
  structure(
    list(table = tab, partition = object),
    class = "summary.icons_partition"
  )
}

#' @export
print.summary.icons_partition <- function(x, ...) {
  print(x$partition)
  if (nrow(x$table)) {
    cat("\n")
    print(x$table, row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' Community membership of a partition
#'
#' Extracts the community label of every variable, in the original variable
#' order.  Unassigned (singleton) variables get `0`, or `NA` if you prefer to
#' treat them as missing.
#'
#' @param x An `"icons_partition"` object from [icons_detect()].
#' @param singleton_na Logical; label singletons `NA` instead of `0`.
#'
#' @return An integer vector of length `p`, named with the variable names when
#'   they are available.
#'
#' @examples
#' data(sim)
#' part <- icons_detect(cor(sim), threshold = 0.6)
#' table(as_membership(part))
#'
#' @export
as_membership <- function(x, singleton_na = FALSE) {
  if (!inherits(x, "icons_partition")) {
    abort_icons("`x` must be an <icons_partition>, not ", class(x)[1L], ".",
                class = "icons_type_error")
  }
  m <- x$membership
  if (singleton_na) m[m == 0L] <- NA_integer_
  names(m) <- x$labels
  m
}

#' Variable indices belonging to given communities
#'
#' @param x An `"icons_partition"` object.
#' @param which Integer vector of community numbers.  Defaults to all of them.
#'
#' @return A list with `index`, the original variable indices in block order,
#'   and `community`, the matching community label.
#'
#' @examples
#' data(sim)
#' part <- icons_detect(cor(sim), threshold = 0.6)
#' block_index(part, which = 1:2)$index
#'
#' @export
block_index <- function(x, which = seq_along(x$sizes)) {
  if (!inherits(x, "icons_partition")) {
    abort_icons("`x` must be an <icons_partition>, not ", class(x)[1L], ".",
                class = "icons_type_error")
  }
  which <- as.integer(which)
  k <- length(x$sizes)
  if (any(which < 1L | which > k)) {
    abort_icons("`which` must be between 1 and ", k, ".",
                class = "icons_value_error")
  }
  end <- cumsum(x$sizes)
  start <- end - x$sizes + 1L
  idx <- unlist(lapply(which, function(j) x$order[start[j]:end[j]]),
                use.names = FALSE)
  list(index = idx, community = rep.int(which, x$sizes[which]))
}

#' Reorder a matrix by a detected partition
#'
#' Permutes rows and columns so community members are adjacent, which is what
#' makes the block structure visible in a heat map.
#'
#' @param W A square matrix with the same dimension as the matrix that produced
#'   `partition`.
#' @param partition An `"icons_partition"` object.
#' @param communities_only Logical; drop the singleton set.
#'
#' @return The reordered matrix.
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#' part <- icons_detect(W, threshold = 0.6)
#' plot_matrix(reorder_matrix(W, part), partition = part)
#'
#' @export
reorder_matrix <- function(W, partition, communities_only = FALSE) {
  if (!inherits(partition, "icons_partition")) {
    abort_icons("`partition` must be an <icons_partition>, not ",
                class(partition)[1L], ".", class = "icons_type_error")
  }
  W <- as_weight_matrix(W, arg = "W", symmetric = FALSE)
  if (ncol(W) != partition$p) {
    abort_icons("`W` has ", ncol(W), " columns but `partition` describes ",
                partition$p, " variables.", class = "icons_dim_error")
  }
  idx <- partition$order
  if (communities_only) idx <- idx[seq_len(sum(partition$sizes))]
  W[idx, idx, drop = FALSE]
}
