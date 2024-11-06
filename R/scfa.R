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
  K <- CID_len - 1
  n <- nrow(data)
  idx_out <- get_index(CID_len, CID, Clist)[[1]]
  CID_temp <- CID[-CID_len]
  clist_in <- get_index(1:K, CID, Clist)[[1]]
  clist_out <- setdiff(1:ncol(data), clist_in)
  Clist_temp <- c(clist_in, clist_out)
  p <- ncol(data) - length(idx_out)
  L <- matrix(0, p, K)

  indices <- rep(1:K, CID_temp)
  L[cbind(1:p, indices)] <- 1

  Y_Data <- data[, clist_in]
  YT_Data <- t(Y_Data)
  YT_mean <- apply(YT_Data, 1, mean)
  YT_Data_Cent <- YT_Data - matrix(rep(YT_mean, n), p, n, byrow = FALSE)

  F_HAT_UB <- solve(t(L) %*% L) %*% (t(L) %*% YT_Data_Cent)

  return(list(loading = L,
              factorscore = F_HAT_UB,
              x = Y_Data))
}
