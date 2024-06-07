#' Parameter Tuning for SigmaU
#'
#' This function tunes parameters for the SigmaU function by evaluating different
#' combinations of lambda and cut-off values.
#'
#' @param Wp A symmetric matrix representing pairwise potentials (weights).
#' @param data A data matrix used in the SigmaU function.
#' @param prctile_vec A vector of percentiles for the pairwise potentials to consider.
#' @param lam_vec A vector of lambda values to test.
#'
#' @return A list containing the optimal lambda value (`lambda_out`) and the corresponding cut-off value (`cut_out`).
#' @export
#'
#' @examples
#' Wp <- matrix(runif(100), nrow = 10)
#' data <- matrix(runif(200), nrow = 20)
#' prctile_vec <- c(10, 20, 30)
#' lam_vec <- c(0.1, 0.5, 1.0)
#' result <- param_tuning_sigmau(Wp, data, prctile_vec, lam_vec)
#' print(result)
#'
param_tuning_sigmau <- function(Wp, data, prctile_vec, lam_vec) {
  # Convert pairwise potential matrix to a vector of non-redundant entries
  nlogp <- get_vectorform(as.matrix(Wp))

  # Calculate cut-off values based on the given percentiles
  gr <- quantile(nlogp, prctile_vec / 100) # percentiles should be in [0, 1]

  nodeLen <- nrow(Wp)
  edgeLen <- length(nlogp)
  sigmau_vec <- matrix(0, nrow = length(lam_vec), ncol = length(gr))

  # Selection of lambda0 begins here
  for (i in seq_along(lam_vec)) {
    lambda <- lam_vec[i]

    for (j in seq_along(gr)) {
      r <- gr[j]

      # Perform greedy clustering with the current cut-off and lambda
      result <- greedy(Wp, r, lambda)
      Clist_temp <- result$Clist
      CID_temp <- result$CID

      # Calculate SigmaU for the current clustering
      sigmau_result <- get_sigmau(data, CID_temp, Clist_temp)
      sigmau_norm <- sigmau_result$sigma_u_norm

      # Store the SigmaU value with a penalty for the number of clusters
      sigmau_vec[i, j] <- sigmau_norm + log(length(CID_temp))
    }
  }

  # Choose the optimal lambda and cut-off
  min_index <- which.min(sigmau_vec)

  # Convert linear index to row and column indices
  row_idx <- ((min_index - 1) %% nrow(sigmau_vec)) + 1
  col_idx <- ((min_index - 1) %/% nrow(sigmau_vec)) + 1

  cut_out <- gr[col_idx]
  lambda_out <- lam_vec[row_idx]

  return(list(lambda_out = lambda_out, cut_out = cut_out))
}
