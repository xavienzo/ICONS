#' Parameter Tuning for Optimal Subnetwork Structures
#'
#' This function tunes parameters to identify optimal subnetwork structures
#' by minimizing the Frobenius norm of the residual covariance matrix.
#' It uses a grid search over different combinations of lambda and cut-off values.
#'
#' @param Wp A symmetric matrix representing pairwise potentials (weights).
#' @param data A data matrix used in the SigmaU function.
#' @param prctile_vec A vector of percentiles for the pairwise potentials to consider.
#' @param lam_vec A vector of lambda values to test.
#' @param method A string indicating the method to use for SigmaU calculation (default is "Sample").
#' @param ncores The number of cores to use for parallel processing (default is NULL, which uses all available cores - 1).
#' @param use_parallel A logical value indicating whether to use parallel processing. Defaults to TRUE.
#'
#' @return A list containing the optimal lambda value (`lambda_out`) and the corresponding cut-off value (`cut_out`).
#' @export
#' @importFrom parallel makeCluster detectCores clusterExport parLapply stopCluster mclapply
#' @importFrom stats quantile

param_tuning_sigmau <- function(Wp, data, prctile_vec, lam_vec, method = "Sample", ncores = NULL, use_parallel = TRUE) {
  # Convert pairwise potential matrix to a vector of non-redundant entries
  nlogp <- get_vectorform(as.matrix(Wp))

  # Calculate cut-off values based on the given percentiles
  gr <- quantile(nlogp, prctile_vec / 100)  # percentiles should be in [0, 1]

  nodeLen <- nrow(Wp)
  edgeLen <- length(nlogp)
  sigmau_vec <- matrix(0, nrow = length(lam_vec), ncol = length(gr))

  # Set up parallel processing (if ncores is provided)
  if (is.null(ncores)) {
    ncores <- detectCores() - 1  # Use one less than the total number of cores
  }

  # Define the function for parameter tuning for a given combination of lambda and cut-off
  tune_param_set <- function(lambda, r) {
    # Perform greedy clustering with the current cut-off and lambda
    result <- dense(Wp, r, lambda)
    Clist_temp <- result$Clist
    CID_temp <- result$CID

    # Calculate SigmaU for the current clustering
    sigmau_result <- scfa(data, CID_temp, Clist_temp, method = method)
    sigmau_norm <- sigmau_result$sigma_u_norm
    sigmau_norm_d <- sigmau_result$sigma_u_diag_norm
    # Store the SigmaU value with a penalty for the number of clusters
    return(sigmau_norm + sigmau_norm_d)
  }

  # Check if parallel processing is enabled
  if (use_parallel) {
    # Check the operating system and set up parallel computing
    if (.Platform$OS.type == "windows") {
      # Windows system: use makeCluster
      cl <- makeCluster(ncores)

      # Export necessary variables to the cluster (if needed)
      clusterExport(cl, list("Wp", "data", "gr", "lam_vec", "dense", "scfa", "method", "tune_param_set"))

      # Parallel processing using parLapply
      results <- parLapply(cl, seq_along(lam_vec), function(i) {
        lambda <- lam_vec[i]
        sapply(gr, function(r) tune_param_set(lambda, r))
      })

      # Stop the cluster after parallel computation
      stopCluster(cl)
    } else {
      # macOS/Linux system: use mclapply (supports multicore)
      results <- mclapply(seq_along(lam_vec), function(i) {
        lambda <- lam_vec[i]
        sapply(gr, function(r) tune_param_set(lambda, r))
      }, mc.cores = ncores)
    }

    # Convert the list of results into a matrix
    sigmau_vec <- do.call(rbind, results)
  } else {
    # Run the tuning process without parallel processing (sequential mode)
    for (i in seq_along(lam_vec)) {
      lambda <- lam_vec[i]
      for (j in seq_along(gr)) {
        r <- gr[j]

        # Perform greedy clustering with the current cut-off and lambda
        result <- dense(Wp, r, lambda)
        Clist_temp <- result$Clist
        CID_temp <- result$CID

        # Calculate SigmaU for the current clustering
        sigmau_result <- scfa(data, CID_temp, Clist_temp, method = method)
        sigmau_norm <- sigmau_result$sigma_u_norm
        sigmau_norm_d <- sigmau_result$sigma_u_diag_norm
        # Store the SigmaU value with a penalty for the number of clusters
        sigmau_vec[i, j] <- sigmau_norm + sigmau_norm_d
      }
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
