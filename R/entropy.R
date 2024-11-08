#' Calculate \eqn{\Sigma_u}{Sigma_u}, the unexplained residual
#'
#' This function calculates \eqn{\Sigma_u}{Sigma_u} and its Frobenius norm,
#' given the data and subnetwork detection results.
#'
#' @param data A data matrix.
#' @param CID A vector of cluster sizes.
#' @param Clist A vector of cluster assignments.
#'
#' @return A list containing \eqn{\Sigma_u}{Sigma_u} (`sigma_u`) and its Frobenius norm (`sigma_u_norm`).
#' @importFrom stats cov quantile
#' @export

get_sigmau <- function(data, CID, Clist) {

  factor_analysis <- scfa(data, CID, Clist)
  F_HAT <- factor_analysis$factorscore
  L <- factor_analysis$loading

  n <- factor_analysis$n.obs
  k <- factor_analysis$n.network
  k0 <- k + 1
  p <- factor_analysis$n.var
  p0 <- ncol(data)

  FT_HAT <- t(F_HAT)

  LT_FULL <- matrix(0, k, p0)
  LT_FULL[, 1:p] <- t(L)

  X <- data[, Clist]
  U <- X - FT_HAT %*% LT_FULL

  SIGMA_U <- cov(U)
  diag(SIGMA_U) <- 0

  SIGMA_UF <- norm(SIGMA_U, type = "F")

  ######################
  if (k > 1) {
    SIGMA_F0 <- diag(diag(cov(FT_HAT)))
  } else {
    SIGMA_F0 <- cov(FT_HAT)
  }
  SIGMA_X <- cov(X)
  L_SIGMAF_LT <- L %*% SIGMA_F0 %*% t(L)

  L_SIGMAF_LT_full <- matrix(0, p0, p0)
  L_SIGMAF_LT_full[1:p, 1:p] <- L_SIGMAF_LT

  SIGMA_U_DIAG <- SIGMA_X - L_SIGMAF_LT_full
  diag(SIGMA_U_DIAG) <- 0
  SIGMA_UF_DIAG <- norm(SIGMA_U_DIAG, type = "F")

  return(list(sigma_u = SIGMA_U,
              sigma_u_norm = SIGMA_UF,
              sigma_u_diag = SIGMA_U_DIAG,
              sigma_u_diag_norm = SIGMA_UF_DIAG))
}
