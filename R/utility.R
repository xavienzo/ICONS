#' Convert a Square-Form Matrix to a Vector
#'
#' Given a distance matrix, this function converts it into a half-vector form,
#' (a vector of the upper triangular part excluding the diagonal).
#'
#' @param dist A square symmetric matrix.
#' @return A vector containing the upper triangular part of the matrix.
#' @export
#' @examples
#' data(sim)
#' get_vectorform(cor(sim))

get_vectorform <- function(dist) {
  # check matrix form
  if (!is.matrix(dist)) {
    stop("Input must be a matrix.")
  }
  if (!identical(dist, t(dist))) {
    stop("Input must be a square symmetric matrix.")
  }
  # return vector form
  return(dist[upper.tri(dist)])
}

#' Calculate \eqn{\Sigma_u}{Sigma_u}, the unexplained residual
#'
#' This function calculates \eqn{\Sigma_u}{Sigma_u} and its Frobenius norm,
#' given the data and subnetwork detection results.
#'
#' @param data A data matrix.
#' @param cid A vector of cluster sizes.
#' @param clist A vector of cluster assignments.
#'
#' @return A list containing the SigmaU matrix (`sigma_u`) and its Frobenius norm (`sigma_u_norm`).
#' @export

get_sigmau <- function(data, cid, clist) {
  K <- length(cid)
  n <- nrow(data)
  p <- ncol(data)
  L <- matrix(0, p, K)

  indices <- rep(1:K, cid)
  L[cbind(1:p, indices)] <- 1

  Y_Data <- data[, clist]
  YT_Data <- t(Y_Data)
  YT_mean <- apply(YT_Data, 1, mean)
  YT_Data_Cent <- YT_Data - matrix(rep(YT_mean, n), p, n, byrow = FALSE)

  F_HAT_UB <- solve(t(L) %*% L) %*% (t(L) %*% YT_Data_Cent)

  sigma_f <- cov(t(F_HAT_UB))
  sigma <- cov(Y_Data)

  lsigmaflt <- L %*% sigma_f %*% t(L)
  sigma_u <- sigma - lsigmaflt

  sigma_u_norm <- norm(sigma_u, type = "F")

  return(list(sigma_u = sigma_u, sigma_u_norm = sigma_u_norm))
}
