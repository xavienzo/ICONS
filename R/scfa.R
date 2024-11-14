#' Semi-Confirmatory Factor Analysis
#'
#' This function performs a semi-confirmatory factor analysis (SCFA) on the provided dataset.
#'
#' @param data A data frame or matrix where rows represent observations and columns represent variables.
#' @param CID A numeric vector indicating the number of variables associated with each factor.
#' @param Clist A numeric vector indicating the column indices of the variables to be used in the analysis.
#' @param method Character. Method of residual term estimation: "Sample" or "MLE".
#'
#' @return A list containing:
#' \describe{
#'   \item{loading}{A matrix of factor loadings.}
#'   \item{factorscore}{A matrix of estimated factor scores.}
#'   \item{sigma_u}{Matrix of off-diagonal residuals (Sample or MLE).}
#'   \item{sigma_u_norm}{Frobenius norm of sigma_u.}
#'   \item{sigma_u_diag}{Diagonalized residuals (Sample or MLE).}
#'   \item{sigma_u_diag_norm}{Frobenius norm of sigma_u_diag.}
#'   \item{n.network}{Integer. Number of factors.}
#'   \item{n.var}{Integer. Number of variables.}
#'   \item{n.obs}{Integer. Number of observations.}
#'   \item{method}{The chosen estimation method.}
#' }
#'
#' @export
#' @importFrom stats cov

scfa <- function(data, CID, Clist, method = "Sample") {

  method <- match.arg(method, choices = c("Sample", "MLE"))

  # parameters
  k0 <- length(CID)
  k <- k0 - 1
  n <- nrow(data)
  CID_temp <- CID[-k0]
  p <- sum(CID_temp)
  p0 <- ncol(data)

  # Set up factor loading matrix
  L <- matrix(0, p, k)
  indices <- rep(1:k, CID_temp)
  L[cbind(1:p, indices)] <- 1
  LT <- t(L)

  # Data centering for selected columns
  X0 <- data[, Clist, drop = FALSE]
  XT0_cent <- t(sweep(X0, 2, colMeans(X0)))
  XT_cent <- XT0_cent[1:p, , drop = FALSE]

  F_HAT <- solve(LT %*% L) %*% (LT %*% XT_cent)

  SIGMA_F <- cov(t(F_HAT))
  L_SIGMAF_LT <- L %*% SIGMA_F %*% LT

  if (method == "Sample") {
    U <- t(XT_cent) - t(F_HAT) %*% LT
    SIGMA_U <- diag(diag(cov(U)))
  } else {
    SIGMA_U <- diag(diag(cov(t(XT_cent)) - L_SIGMAF_LT))
  }

  F_HAT_FINAL <- solve(LT %*% solve(SIGMA_U) %*% L) %*% LT %*% solve(SIGMA_U) %*% XT_cent

  ################ GET SIGMAU #################
  SIGMA_X0 <- cov(t(XT0_cent))
  SIGMA_F <- cov(t(F_HAT_FINAL))
  L_SIGMAF_LT_full <- matrix(0, p0, p0)
  L_SIGMAF_LT_full[1:p, 1:p] <- L %*% SIGMA_F %*% t(L)

  SIGMA_U_OFF <- SIGMA_X0 - L_SIGMAF_LT_full
  diag(SIGMA_U_OFF) <- 0
  SIGMA_UF_OFF <- norm(SIGMA_U_OFF, type = "F")

  ################ GET SIGMAU 0 (DIAGNAL) #################
  SIGMA_F0 <- if (k > 1) diag(diag(SIGMA_F)) else SIGMA_F
  L_SIGMAF_LT_full[1:p, 1:p] <- L %*% SIGMA_F0 %*% t(L)

  SIGMA_U_DIAG <- SIGMA_X0 - L_SIGMAF_LT_full
  diag(SIGMA_U_DIAG) <- 0
  SIGMA_UF_DIAG <- norm(SIGMA_U_DIAG, type = "F")

  return(list(loading = L,
              factorscore = F_HAT,
              sigma_u = SIGMA_U_OFF,
              sigma_u_norm = SIGMA_UF_OFF,
              sigma_u_diag = SIGMA_U_DIAG,
              sigma_u_diag_norm = SIGMA_UF_DIAG,
              n.network = k,
              n.var = p,
              n.obs = n,
              method = method))
}
