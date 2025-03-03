#' Semi-Confirmatory Factor Analysis
#'
#' This function performs a semi-confirmatory factor analysis (SCFA) on the provided dataset.
#'
#' @param data A scaled data frame or matrix where rows represent observations and columns represent variables.
#' @param CID A numeric vector indicating the number of variables associated with each factor.
#' @param Clist A numeric vector indicating the column indices of the variables to be used in the analysis.
#' @param method Character. Method of residual term estimation: "Sample" or "MLE".
#' @param epsilon A small regularization value for solving the inverse of Sigma_U.
#' @param remove.singletons Logical. Remove the singleton set for SCFA or not.
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

scfa <- function(data,
                  CID,
                  Clist,
                  method = "Sample",
                  epsilon = 1e-6,
                  remove.singletons = TRUE) {

  method <- match.arg(method, choices = c("Sample", "MLE"))

  # parameters
  if (remove.singletons == T) {
    k <- length(CID) - 1
    CID_temp <- CID[-length(CID)]
  } else {
    k <- length(CID)
    CID_temp <- CID
  }
  n <- nrow(data)
  p <- sum(CID_temp)
  p0 <- ncol(data)

  # Set up factor loading matrix
  L <- matrix(0, p, k)
  indices <- rep(1:k, CID_temp)
  L[cbind(1:p, indices)] <- 1
  LT <- t(L)

  # Data centering for selected columns
  X0 <- data[, Clist, drop = FALSE]
  XT <- t(X0)[1:p, , drop = FALSE]

  if (length(CID_temp) > 1) {
    F_HAT <- diag(1 / CID_temp) %*% (LT %*% XT)
  } else {
    F_HAT <- (1 / CID_temp) * (LT %*% XT)  # Use `*` for scalar multiplication
  }

  SIGMA_F <- cov(t(F_HAT))
  L_SIGMAF_LT <- L %*% SIGMA_F %*% LT

  if (method == "Sample") {
    U <- t(XT) - t(F_HAT) %*% LT
    SIGMA_U <- diag(diag(cov(U)))
  } else {
    SIGMA_U <- diag(diag(cov(t(XT)) - L_SIGMAF_LT))
  }

  INV_SIGMA_U <- tryCatch({
    diag(1/diag(SIGMA_U))
  }, error = function(e) {
    SIGMA_U_reg <- SIGMA_U + diag(epsilon, nrow(SIGMA_U))
    diag(1/diag(SIGMA_U_reg))
  })

  if (length(CID_temp) > 1) {
    F_HAT_FINAL <- diag(1 / tapply(diag(INV_SIGMA_U), rep(seq_along(CID_temp), CID_temp), sum)) %*% LT %*% INV_SIGMA_U %*% XT
  } else {
    F_HAT_FINAL <- solve(LT %*% INV_SIGMA_U %*% L) %*% LT %*% INV_SIGMA_U %*% XT
  }

  ################ GET SIGMAU #################
  SIGMA_F <- cov(t(F_HAT_FINAL))
  LT_FULL <- matrix(0, k, p0)
  LT_FULL[, 1:p] <- t(L)
  U <- X0 - t(F_HAT_FINAL) %*% LT_FULL
  SIGMA_U_OFF <- cov(U)
  diag(SIGMA_U_OFF) <- 0
  SIGMA_UF_OFF <- norm(SIGMA_U_OFF, type = "F")

  ################ GET SIGMAU 0 (DIAGNAL) #################
  SIGMA_F0 <- if (k > 1) diag(diag(SIGMA_F)) else SIGMA_F
  L_SIGMAF_LT_full <- matrix(0, p0, p0)
  L_SIGMAF_LT_full[1:p, 1:p] <- L %*% SIGMA_F0 %*% t(L)

  SIGMA_U_DIAG <- cov(X0) - L_SIGMAF_LT_full
  diag(SIGMA_U_DIAG) <- 0
  SIGMA_UF_DIAG <- norm(SIGMA_U_DIAG, type = "F")

  return(list(loading = L,
              factorscore = F_HAT_FINAL,
              sigma_u = SIGMA_U_OFF,
              sigma_u_norm = SIGMA_UF_OFF,
              sigma_u_diag = SIGMA_U_DIAG,
              sigma_u_diag_norm = SIGMA_UF_DIAG,
              n.network = k,
              n.var = p,
              n.obs = n,
              method = method))
}
