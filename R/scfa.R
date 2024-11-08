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

scfa <- function(data, CID, Clist) {
  CID_len <- length(CID)
  k <- CID_len - 1
  n <- nrow(data)

  CID_temp <- CID[-CID_len]
  clist_in <- get_index(1:k, CID, Clist)[[1]]
  p <- length(clist_in)

  L <- matrix(0, p, k)
  indices <- rep(1:k, CID_temp)
  L[cbind(1:p, indices)] <- 1

  X <- data[, clist_in]
  XT <- t(X)
  XT_mean <- apply(XT, 1, mean)
  XT_cent <- XT - matrix(rep(XT_mean, n), p, n, byrow = FALSE)

  F_HAT <- solve(t(L) %*% L) %*% (t(L) %*% XT_cent)

  return(list(loading = L,
              factorscore = F_HAT,
              x = X,
              n.network = k,
              n.var = p,
              n.obs = n))
}
