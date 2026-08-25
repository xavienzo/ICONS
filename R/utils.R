## Internal helpers -----------------------------------------------------------
## Argument checking is centralised here so that every entry point reports
## problems the same way and names the offending argument.

stop_icons <- function(..., class = NULL, call = sys.call(-1L)) {
  msg <- paste0(..., collapse = "")
  structure(
    class = c(class, "icons_error", "error", "condition"),
    list(message = msg, call = call)
  )
}

abort_icons <- function(..., class = NULL) {
  stop(stop_icons(..., class = class, call = sys.call(-1L)))
}

# Coerce to a plain numeric matrix and check the properties every entry point
# needs.  `symmetric = TRUE` additionally requires (approximate) symmetry.
as_weight_matrix <- function(W, arg = "W", symmetric = TRUE, tol = 1e-8) {
  if (is.data.frame(W)) W <- as.matrix(W)
  if (!is.matrix(W)) {
    abort_icons("`", arg, "` must be a matrix or data frame, not ",
                class(W)[1L], ".", class = "icons_type_error")
  }
  if (!is.numeric(W)) {
    abort_icons("`", arg, "` must be numeric.", class = "icons_type_error")
  }
  if (nrow(W) != ncol(W)) {
    abort_icons("`", arg, "` must be square; it is ", nrow(W), " by ",
                ncol(W), ".", class = "icons_dim_error")
  }
  if (nrow(W) < 2L) {
    abort_icons("`", arg, "` must have at least 2 rows.",
                class = "icons_dim_error")
  }
  if (anyNA(W)) {
    abort_icons("`", arg, "` must not contain missing values.",
                class = "icons_na_error")
  }
  if (!all(is.finite(W))) {
    abort_icons("`", arg, "` must contain only finite values.",
                class = "icons_na_error")
  }
  if (symmetric) {
    d <- max(abs(W - t(W)))
    if (d > tol * max(1, max(abs(W)))) {
      abort_icons("`", arg, "` must be symmetric; the largest asymmetry is ",
                  format(d, digits = 3), ".", class = "icons_sym_error")
    }
    # Enforce exact symmetry so downstream results do not depend on rounding.
    W <- (W + t(W)) / 2
  }
  storage.mode(W) <- "double"
  W
}

as_data_matrix <- function(data, arg = "data") {
  if (is.data.frame(data)) {
    num <- vapply(data, is.numeric, logical(1L))
    if (!all(num)) {
      abort_icons("`", arg, "` must be all-numeric; column(s) ",
                  paste(names(data)[!num], collapse = ", "),
                  " are not.", class = "icons_type_error")
    }
    data <- as.matrix(data)
  }
  if (!is.matrix(data) || !is.numeric(data)) {
    abort_icons("`", arg, "` must be a numeric matrix or data frame.",
                class = "icons_type_error")
  }
  if (anyNA(data)) {
    abort_icons("`", arg, "` must not contain missing values.",
                class = "icons_na_error")
  }
  storage.mode(data) <- "double"
  data
}

check_scalar_number <- function(x, arg, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    abort_icons("`", arg, "` must be a single finite number.",
                class = "icons_type_error")
  }
  if (x < lower || x > upper) {
    abort_icons("`", arg, "` must be between ", lower, " and ", upper,
                "; got ", x, ".", class = "icons_value_error")
  }
  as.numeric(x)
}

check_count <- function(x, arg, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x)) {
    abort_icons("`", arg, "` must be a single whole number.",
                class = "icons_type_error")
  }
  if (x < lower) {
    abort_icons("`", arg, "` must be at least ", lower, "; got ", x, ".",
                class = "icons_value_error")
  }
  as.integer(x)
}

# Frobenius norm without materialising a copy of the matrix.
fro <- function(x) sqrt(sum(x * x))

#' Half-vectorisation of a symmetric matrix
#'
#' Returns the strict upper triangle of a symmetric matrix as a vector, the
#' form used for thresholding edge weights by quantile.
#'
#' @param x A square symmetric numeric matrix.
#' @param check Logical; verify symmetry. Set to `FALSE` to skip the check on
#'   very large matrices.
#'
#' @return A numeric vector of length `ncol(x) * (ncol(x) - 1) / 2`.
#'
#' @examples
#' data(sim)
#' v <- half_vec(cor(sim))
#' quantile(v, c(0.95, 0.99))
#'
#' @export
half_vec <- function(x, check = TRUE) {
  x <- as_weight_matrix(x, arg = "x", symmetric = check)
  x[upper.tri(x)]
}
