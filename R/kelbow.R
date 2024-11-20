#' Determine the Elbow Point for k in SCFA
#'
#' This function computes the residual error for different numbers of factors (k) in a semi-confirmatory
#' factor analysis (SCFA), helping to identify the optimal number of factors using the "elbow" method.
#'
#' @param data A data frame or matrix where rows represent observations and columns represent variables.
#' @param CID A numeric vector indicating the number of variables associated with each factor.
#' @param Clist A numeric vector indicating the column indices of the variables reordered by factors.
#' @param k Integer. The maximum number of factors to evaluate. Defaults to the length of `CID`.
#' @param method Character. The method of residual term estimation: "Sample" or "MLE". Defaults to "Sample".
#' @param epsilon Numeric. A small regularization term added to diagonal elements for numerical stability. Defaults to 1e-6.
#'
#' @return A numeric vector where each element represents the residual norm (sum of off-diagonal and diagonal norms)
#' for the corresponding number of factors (from 0 to `k`).
#'
#' @export
#'
k.elbow <- function(data,
                    CID,
                    Clist,
                    k = length(CID),
                    method = "Sample",
                    epsilon = 1e-6) {

  # Validate and initialize parameters
  k <- min(k, length(CID)) # Ensure k is within the length of CID
  covmat <- cov(data)
  diag(covmat) <- 0
  fullgraph_norm <- 2 * norm(covmat, type = "F")
  sigmau_k_f_vector <- numeric(k + 1)
  sigmau_k_f_vector[1] <- fullgraph_norm

  # Compute residual norms for increasing numbers of factors
  for (kk in 1:k) {
    clist_in <- Clist[1:sum(CID[1:kk])]
    CID_temp <- c(CID[1:kk], ncol(data) - length(clist_in))
    sigmau <- scfa(data, CID_temp, clist_in, method, epsilon)
    sigmau_k_f_vector[kk + 1] <- sigmau$sigma_u_norm + sigmau$sigma_u_diag_norm
  }

  return(sigmau_k_f_vector)
}
