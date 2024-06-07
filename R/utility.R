#' Convert a Squareform Matrix to a Vector
#'
#' Given a squareform matrix, this function converts it into a vector containing
#' the lower triangular part of the matrix (excluding the diagonal).
#'
#' @param mat A square-form matrix.
#' @return A vector containing the lower triangular part of the matrix.
#' @export
#' @examples
#' data(sim)
#' squareform(sim)
#'
#' mat <- matrix(1:6, nrow = 3)
#' squareform(mat)

squareform <- function(mat) {
  indices <- lower.tri(mat)
  result <- mat[indices]
  return(result)
}

#' Compute SigmaU and its Norm
#'
#' This function calculates the SigmaU matrix and its Frobenius norm for given data and clustering results.
#'
#' @param data A data matrix.
#' @param cid A vector of cluster sizes.
#' @param clist A vector of cluster assignments.
#'
#' @return A list containing the SigmaU matrix (`sigma_u`) and its Frobenius norm (`sigma_u_norm`).
#' @export
#'
#' @examples
#' data <- matrix(runif(200), nrow = 20)
#' cid <- c(5, 5)
#' clist <- sample(1:2, 20, replace = TRUE)
#' result <- sigmau(data, cid, clist)
#' print(result)
#'
sigmau <- function(data, cid, clist) {
  K <- length(cid)
  n <- nrow(data)
  p <- ncol(data)
  L <- matrix(0, p, K)

  indices <- rep(1:K, cid)
  L[cbind(1:p, indices)] <- 1

  Y_Data <- data[, clist]
  YT_Data <- t(Y_Data)
  YT_mean <- apply(YT_Data, 1, mean)
  # YT_Data_Cent <- sweep(YT_Data, 2, YT_mean)
  YT_Data_Cent <- YT_Data - matrix(rep(YT_mean, n), p, n, byrow = FALSE)

  F_HAT_UB <- solve(t(L) %*% L) %*% (t(L) %*% YT_Data_Cent)

  sigma_f <- cov(t(F_HAT_UB))
  sigma <- cov(Y_Data)

  lsigmaflt <- L %*% sigma_f %*% t(L)
  sigma_u <- sigma - lsigmaflt

  sigma_u_norm <- norm(sigma_u, type = "F")

  return(list(sigma_u = sigma_u, sigma_u_norm = sigma_u_norm))
}
