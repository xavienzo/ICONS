#' Semi-Confirmatory Factor Analysis
#'
#' This function performs a semi-confirmatory factor analysis (SCFA) on the provided dataset.
#'
#' @param data A data frame or matrix where rows represent observations and columns represent variables.
#' @param cid A numeric vector indicating the number of variables associated with each factor.
#' @param clist A numeric vector indicating the column indices of the variables to be used in the analysis.
#'
#' @return A list containing:
#' \item{loading}{A matrix of factor loadings.}
#' \item{factorscore}{A matrix of estimated factor scores.}
#'
#' @export

scfa <- function(data, cid, clist) {
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

  return(list(loading = L,
              factorscore = F_HAT_UB,
              x = Y_Data))
}
