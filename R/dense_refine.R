#' Get Indices of Submatrices from a Block-Structured Adjacency Matrix
#'
#' This function extracts the indices of the elements corresponding to one or more submatrices
#' (blocks) in a larger adjacency matrix. The submatrices are defined by a vector containing the
#' number of elements in each block, and a vector that holds the rearranged index of the full matrix
#' based on the submatrix structure.
#'
#' @param blockid A numeric value or vector indicating the block(s) (submatrix/submatrices) to be selected.
#'                Block IDs should be integers corresponding to positions in the `CID` vector.
#' @param CID A numeric vector where each element represents the size (number of elements) of each submatrix.
#' @param Clist A numeric vector indicating the rearranged index of the full adjacency matrix based on the
#'              submatrices. The length of `Clist` should be equal to the sum of the values in `CID`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{indices}{A numeric vector of the indices in `Clist` corresponding to the specified block(s).}
#'   \item{block_ids}{A numeric vector containing the block IDs for each corresponding index in `indices`.}
#' }
#'
#' @examples
#' CID <- c(3, 4, 2)  # 3 nodes in block 1, 4 nodes in block 2, 2 nodes in block 3
#' Clist <- c(3, 1, 2, 6, 4, 5, 7, 9, 8)  # Rearranged index of the full matrix
#'
#' # Get indices for block 1 and 2
#' result <- get_index(blockid = c(1, 2), CID, Clist)
#' print(result$indices)    # Indices of the rearranged matrix
#' print(result$block_ids)  # Corresponding block ids
#'
#' @export

get_index <- function(blockid, CID, Clist){
  # safeguard for input like 1:2
  blockid <- as.numeric(blockid)
  # Calculate start and end positions for each block based on CID
  start_pos <- cumsum(c(0, head(CID, -1))) + 1 # Start positions for each block
  end_pos <- cumsum(CID) # End positions for each block

  # Create a vector of logicals to identify valid blockids
  valid_blockids <- blockid[blockid <= length(CID) & blockid >= 1]

  if (length(valid_blockids) != length(blockid)) {
    stop("One or more block IDs are out of range.")
  }

  # Get the start and end indices for each block in blockid
  block_indices <- mapply(function(start, end) Clist[start:end],
                          start_pos[valid_blockids], end_pos[valid_blockids])

  # Flatten the list of indices into a single vector
  final_indices <- as.vector(unlist(block_indices))

  # Repeat the block ID for the length of each block's indices
  final_block_ids <- rep(valid_blockids, times = CID[valid_blockids])

  return(list(indices = final_indices, block_ids = final_block_ids))
}

dense_refine <- function(W_original, blockid, CID, Clist,
                         threshold = 0.5, lambda = 0.5,
                         param_tuning = F, ...){
  # Index of the original matrix that needs refinement
  idx_refine <- get_index(blockid, CID, Clist)[[1]]
  idx_preserve <- get_index(setdiff(1:length(CID), blockid), CID, Clist)[[1]]

  # Matrix that needs refinement
  W_refine <- W_original[idx_refine, idx_refine]

  # Apply dense algorithm on the new matrix
  if (param_tuning) {   # Apply parameter tuning function if indicated
    tuning_params <- list(...)
    working_params <- param_tuning_sigmau(W_refine, data, prctile_vec, lam_vec)
    threshold <- working_params$cut_out
    lambda <- working_params$lambda_out
  }
  working_results <- dense(W_refine, threshold, lambda)

  # Dense index order by original matrix index
  W_dense_refine <- working_results$W_dense
  idx_refine_reorder <- idx_refine[working_results$Clist]

  # Insert the refined dense submatrix back into the original matrix
  W_dense_full <- W_original[Clist, Clist] -
  idx_replace_start <- ifelse(blockid[1] == 1, 1, as.numeric(CID[blockid[1] - 1]) + 1)
  idx_replace_end <- sum(CID[1:blockid[length(blockid)]])
  idx_replace <- idx_replace_start : idx_replace_end
  W_dense_full[idx_replace, idx_replace] <- W_dense_refine

  if (blockid[1] == 1) {
    idx_full <- c(idx_refine_reorder, idx_preserve)
  } else if (blockid[length(blockid)] != length(CID)) {
    # Case where blockid is not at the start or end
    idx_full <- c(
      Clist[1:cumsum(CID)[blockid[1] - 1]],
      idx_refine_reorder,
      Clist[(cumsum(c(0, head(CID, -1)))[blockid[length(blockid)] + 1] + 1) : length(Clist)]
    )
  } else {
    # Case where blockid ends with the last block
    idx_full <- c(Clist[1:cumsum(CID)[blockid[1] - 1]], idx_refine_reorder)
  }

  return(list(
    W_dense_refine = W_dense_refine,
    W_dense_full = W_dense_full,
    idx_refine = idx_refine,
    Clist_refine = idx_refine_reorder,
    Clist_full = idx_full,
    CID_refine = working_results$CID#,
    #CID_full =
  ))
}
