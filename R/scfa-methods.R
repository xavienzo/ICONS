## S3 methods for "scfa" ------------------------------------------------------

#' @export
print.scfa <- function(x, digits = 4L, ...) {
  cat("<Semi-confirmatory factor analysis>\n")
  cat("  observations :", x$n, "\n")
  cat("  variables    : ", x$p, " (", x$p_modelled, " modelled, ",
      x$n_singletons, if (x$n_singletons == 1L) " singleton)" else " singletons)",
      "\n", sep = "")
  cat("  factors      :", x$K, "\n")
  cat("  Sigma_u      :", x$sigma_u, "\n")
  if (!is.na(x$criterion$value)) {
    cat("  ", x$criterion$criterion, " criterion: ",
        format(x$criterion$value, digits = digits),
        "  (relative ", format(x$criterion$relative, digits = digits), ")\n",
        sep = "")
  }
  invisible(x)
}

#' @export
summary.scfa <- function(object, level = 0.95, ...) {
  ci <- confint(object, level = level)
  structure(
    list(fit = object, coefficients = ci, level = level),
    class = "summary.scfa"
  )
}

#' @export
print.summary.scfa <- function(x, digits = 4L, ...) {
  print(x$fit, digits = digits)
  cat("\nCommunity sizes:\n")
  print(stats::setNames(x$fit$sizes, names(x$fit$a)))
  cat("\nFactor covariance (Sigma_f):\n")
  print(round(x$fit$Sigma_f, digits))
  cat("\nParameter estimates with ", format(100 * x$level), "% Wald intervals:\n",
      sep = "")
  print(round(x$coefficients, digits))
  cat("\nExact standard errors from Theorem 3 of Yang et al. (2024).\n")
  invisible(x)
}

#' Extract SCFA parameter estimates
#'
#' @param object An `"scfa"` object.
#' @param type Which parameters to return: `"all"` stacks the error variances
#'   \eqn{a_{kk}} and the unique factor covariances \eqn{b_{kk'}} into a named
#'   vector, `"a"` returns the error variances, `"Sigma_f"` the factor
#'   covariance matrix.
#' @param ... Ignored.
#'
#' @return A named numeric vector, or a matrix for `type = "Sigma_f"`.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' head(coef(fit))
#'
#' @export
coef.scfa <- function(object, type = c("all", "a", "Sigma_f"), ...) {
  type <- match.arg(type)
  switch(type,
    a = object$a,
    Sigma_f = object$Sigma_f,
    all = {
      pn <- param_names(object$K)
      stats::setNames(c(object$a, object$Sigma_f[upper.tri(object$Sigma_f, TRUE)]),
                      pn)
    }
  )
}

param_names <- function(K) {
  fk <- paste0("F", seq_len(K))
  ij <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  c(paste0("a[", fk, "]"),
    paste0("b[", fk[ij[, "row"]], ",", fk[ij[, "col"]], "]"))
}

#' Variance-covariance of the SCFA parameter estimates
#'
#' The exact variances of Theorem 3 in Yang, Ma, Bi and Chen (2024).  The
#' estimators of distinct blocks are uncorrelated, so the matrix is diagonal.
#'
#' @param object An `"scfa"` object.
#' @param ... Ignored.
#'
#' @return A diagonal matrix over the parameters returned by `coef()`.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' sqrt(diag(vcov(fit)))[1:4]
#'
#' @export
vcov.scfa <- function(object, ...) {
  ut <- upper.tri(object$se$B, diag = TRUE)
  v <- c(object$se$a^2, object$se$B[ut]^2)
  pn <- param_names(object$K)
  structure(diag(v, length(v)), dimnames = list(pn, pn))
}

#' Wald confidence intervals for SCFA parameters
#'
#' @param object An `"scfa"` object.
#' @param parm Optional character or integer subset of the parameters.
#' @param level Confidence level.
#' @param ... Ignored.
#'
#' @return A matrix with columns `estimate`, `se`, and the interval bounds.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' confint(fit, parm = 1:3)
#'
#' @export
confint.scfa <- function(object, parm, level = 0.95, ...) {
  est <- coef(object, "all")
  se <- sqrt(diag(vcov(object)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(estimate = est, se = se,
               lower = est - z * se, upper = est + z * se)
  colnames(out)[3:4] <- paste0(format(100 * c((1 - level) / 2,
                                              1 - (1 - level) / 2),
                                      trim = TRUE, digits = 3), " %")
  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' Factor loading matrix of an SCFA fit
#'
#' Builds \eqn{\hat L = \mathrm{Bdiag}(1_{p_1}, \ldots, 1_{p_K})} in the
#' original variable order.  It is not stored with the fit because the
#' membership determines it exactly, and at large `p` this matrix would dwarf
#' everything else in the object.
#'
#' Named `factor_loadings()` rather than `loadings()` because `stats::loadings`
#' is not a generic --- it is defined as `x$loadings` --- so no S3 method for
#' it is possible, and masking it would be worse than a distinct name.
#'
#' @param object An `"scfa"` object.
#' @param sparse Logical; return a two-column `variable`/`factor` index instead
#'   of the dense matrix.  Useful when `p` is large.
#'
#' @return A `p` by `K` matrix with one non-zero entry per modelled row, or a
#'   data frame when `sparse = TRUE`.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' dim(factor_loadings(fit))
#' head(factor_loadings(fit, sparse = TRUE))
#'
#' @export
factor_loadings <- function(object, sparse = FALSE) {
  if (!inherits(object, "scfa")) {
    abort_icons("`object` must be an <scfa> fit, not ", class(object)[1L], ".",
                class = "icons_type_error")
  }
  idx <- which(object$membership > 0L)
  if (sparse) {
    return(data.frame(variable = idx,
                      name = if (is.null(object$var_names)) NA_character_
                             else object$var_names[idx],
                      factor = object$membership[idx],
                      loading = 1))
  }
  L <- matrix(0, object$p, object$K,
              dimnames = list(object$var_names, names(object$a)))
  L[cbind(idx, object$membership[idx])] <- 1
  L
}

#' @export
nobs.scfa <- function(object, ...) object$n

#' Factor scores for new observations
#'
#' Applies the fitted structure to new data: each score is the mean of the
#' community's variables after subtracting the training centres.
#'
#' @param object An `"scfa"` object.
#' @param newdata A numeric matrix or data frame with the same variables, in
#'   the same order, as the training data.  Defaults to returning the training
#'   scores.
#' @param ... Ignored.
#'
#' @return An `nrow(newdata)` by `K` matrix of factor scores.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' dim(predict(fit, sim[1:5, ]))
#'
#' @export
predict.scfa <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$scores)
  X <- as_data_matrix(newdata, arg = "newdata")
  if (ncol(X) != object$p) {
    abort_icons("`newdata` has ", ncol(X), " columns but the fit used ",
                object$p, ".", class = "icons_dim_error")
  }
  bs <- .block_stats_cpp(X, object$center, object$membership, object$K)
  out <- sweep(bs$T, 2L, object$sizes, "/")
  dimnames(out) <- list(rownames(X), names(object$a))
  out
}

#' Residuals of an SCFA fit
#'
#' The unexplained part \eqn{U = X_c - \hat F \hat L^\top}.  The data are
#' required because the fit does not retain them.
#'
#' @param object An `"scfa"` object.
#' @param data The data the model was fitted to.
#' @param ... Ignored.
#'
#' @return An `n` by `p` matrix.  Singleton columns are returned centred but
#'   otherwise untouched, since no factor was fitted to them.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' U <- residuals(fit, sim)
#' round(range(colMeans(U)), 10)
#'
#' @export
residuals.scfa <- function(object, data, ...) {
  if (missing(data)) {
    abort_icons("`data` is required: an <scfa> fit does not store the data.",
                class = "icons_value_error")
  }
  X <- as_data_matrix(data, arg = "data")
  if (ncol(X) != object$p || nrow(X) != object$n) {
    abort_icons("`data` must be ", object$n, " by ", object$p, ".",
                class = "icons_dim_error")
  }
  U <- sweep(X, 2L, object$center, check.margin = FALSE)
  idx <- object$membership > 0L
  U[, idx] <- U[, idx] - object$scores[, object$membership[idx], drop = FALSE]
  U
}

#' Model-implied covariance matrix
#'
#' Reconstructs \eqn{\hat\Sigma = \hat L\hat\Sigma_f\hat L^\top + \hat\Sigma_u}
#' explicitly.  This is the one function in the package that materialises a
#' \eqn{p \times p} matrix, so it refuses to run for large `p` unless you ask
#' it to.
#'
#' @param object An `"scfa"` object.
#' @param max_p Guard rail: refuse when `p` exceeds this. Defaults to 5000,
#'   about 200 MB.
#' @param ... Ignored.
#'
#' @return A `p` by `p` matrix in the original variable order.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' dim(fitted(fit))
#'
#' @export
fitted.scfa <- function(object, max_p = 5000L, ...) {
  p <- object$p
  if (p > max_p) {
    abort_icons("p = ", p, " would need a ", p, " by ", p,
                " matrix (about ", round(8 * p^2 / 2^20), " MB). ",
                "Raise `max_p` if you really want it.",
                class = "icons_size_error")
  }
  m <- object$membership
  idx <- m > 0L
  Sigma <- matrix(0, p, p, dimnames = list(object$var_names, object$var_names))
  Sigma[idx, idx] <- object$Sigma_f[m[idx], m[idx]]
  du <- numeric(p)
  du[idx] <- object$a[m[idx]]
  if (!is.null(object$Sigma_u_diag)) {
    du <- object$Sigma_u_diag
  } else {
    du[!idx] <- NA_real_
  }
  diag(Sigma) <- diag(Sigma) + du
  Sigma
}

#' Error variances of an SCFA fit
#'
#' The diagonal of \eqn{\hat\Sigma_u}, one entry per variable.  Under the
#' default `sigma_u = "mle"` every variable in community \eqn{k} shares the
#' value \eqn{\hat a_{kk}}.
#'
#' @param object An `"scfa"` object.
#'
#' @return A named numeric vector of length `p`; `NA` for singletons under the
#'   MLE, since no error variance is identified for an unmodelled variable.
#'
#' @examples
#' data(sim)
#' fit <- scfa(sim, icons_detect(cor(sim), threshold = 0.6))
#' head(sigma_u(fit))
#'
#' @export
sigma_u <- function(object) {
  if (!inherits(object, "scfa")) {
    abort_icons("`object` must be an <scfa> fit, not ", class(object)[1L], ".",
                class = "icons_type_error")
  }
  if (!is.null(object$Sigma_u_diag)) return(object$Sigma_u_diag)
  out <- rep(NA_real_, object$p)
  idx <- object$membership > 0L
  out[idx] <- object$a[object$membership[idx]]
  names(out) <- object$var_names
  out
}
