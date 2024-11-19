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

  # Flatten the combinations of lambda and r
  param_combinations <- expand.grid(lambda = lam_vec, r = gr)

  # Define the function for parameter tuning for a given combination of lambda and cut-off
  tune_param_set <- function(lambda, r) {
    result <- dense(Wp, r, lambda)
    sigmau_result <- scfa(data, result$CID, result$Clist, method = method)
    sigmau_norm <- sigmau_result$sigma_u_norm
    sigmau_norm_d <- sigmau_result$sigma_u_diag_norm
    return(list(lambda = lambda, r = r, result = sigmau_norm + sigmau_norm_d))
  }

  # Check if parallel processing is enabled
  if (use_parallel) {
    # Set up parallel processing (if ncores is provided)
    if (is.null(ncores)) {
      ncores <- detectCores() - 1  # Use one less than the total number of cores
    }
    # Check the operating system and set up parallel computing
    if (.Platform$OS.type == "windows") {
      # Windows system: use makeCluster
      cl <- makeCluster(ncores)
      clusterExport(cl, list("Wp", "data", "gr", "lam_vec", "dense", "scfa", "method", "tune_param_set"))

      # Apply the function in parallel
      results <- parLapply(cl, seq_len(nrow(param_combinations)), function(idx) {
        lambda <- param_combinations$lambda[idx]
        r <- param_combinations$r[idx]
        result <- tune_param_set(lambda, r)
        return(result)
      })

      # Stop the cluster
      stopCluster(cl)

      } else {
      # macOS/Linux system: use mclapply (supports multicore)
        results <- mclapply(seq_len(nrow(param_combinations)), function(idx) {
          lambda <- param_combinations$lambda[idx]
          r <- param_combinations$r[idx]
          result <- tune_param_set(lambda, r)
          return(result)
        }, mc.cores = ncores)
      }
  } else {
    # Run the tuning process without parallel processing (sequential mode)
    results <- lapply(seq_len(nrow(param_combinations)), function(idx) {
      lambda <- param_combinations$lambda[idx]
      r <- param_combinations$r[idx]
      result <- tune_param_set(lambda, r)
      return(result)
    })
  }

  # flatten the list into a data frame
  all_results_df <- do.call(rbind, lapply(results, function(res) {
    data.frame(lambda = res$lambda, r = res$r, result = res$result)
  }))
  row.names(all_results_df) <- NULL
  colnames(all_results_df) <- c("Lambda", "CutOff", "SigmaU")

  # optimal parameters
  min_index <- which.min(all_results_df$SigmaU)
  optimal_lambda <- all_results_df$Lambda[min_index]
  optimal_r <- all_results_df$CutOff[min_index]

  return(list(lambda_out = optimal_lambda, cut_out = optimal_r, all = all_results_df))
}
