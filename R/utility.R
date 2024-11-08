#' Convert a Square-Form Matrix to a Vector
#'
#' Given a distance matrix, this function converts it into a half-vector form,
#' (a vector of the upper triangular part excluding the diagonal).
#'
#' @param dist A square symmetric matrix.
#' @return A vector containing the upper triangular part of the matrix.
#' @export
#' @examples
#' data(sim)
#' get_vectorform(cor(sim))

get_vectorform <- function(dist) {
  # check matrix form
  if (!is.matrix(dist)) {
    stop("Input must be a matrix.")
  }
  if (!identical(dist, t(dist))) {
    stop("Input must be a square symmetric matrix.")
  }
  # return vector form
  return(dist[upper.tri(dist)])
}

#' Calculate \eqn{\Sigma_u}{Sigma_u}, the unexplained residual
#'
#' This function calculates \eqn{\Sigma_u}{Sigma_u} and its Frobenius norm,
#' given the data and subnetwork detection results.
#'
#' @param data A data matrix.
#' @param cid A vector of cluster sizes.
#' @param clist A vector of cluster assignments.
#'
#' @return A list containing the SigmaU matrix (`sigma_u`) and its Frobenius norm (`sigma_u_norm`).
#' @importFrom stats cov quantile
#' @export

get_sigmau <- function(data, CID, Clist) {

  factor_analysis <- scfa(data, CID, Clist)
  F_HAT <- factor_analysis$factorscore
  L <- factor_analysis$loading

  n <- factor_analysis$n.obs
  k <- factor_analysis$n.network
  k0 <- k + 1
  p <- factor_analysis$n.var
  p0 <- ncol(data)

  FT_HAT <- t(F_HAT)

  LT_FULL <- matrix(0, k, p0)
  LT_FULL[, 1:p] <- t(L)

  X <- data[, Clist]
  U <- X - FT_HAT %*% LT_FULL

  SIGMA_U <- cov(U)
  diag(SIGMA_U) <- 0

  SIGMA_UF <- norm(SIGMA_U, type = "F")

  ######################
  if (k > 1) {
    SIGMA_F0 <- diag(diag(cov(FT_HAT)))
  } else {
    SIGMA_F0 <- cov(FT_HAT)
  }
  SIGMA_X <- cov(X)
  L_SIGMAF_LT <- L %*% SIGMA_F0 %*% t(L)

  L_SIGMAF_LT_full <- matrix(0, p0, p0)
  L_SIGMAF_LT_full[1:p, 1:p] <- L_SIGMAF_LT

  SIGMA_U_DIAG <- SIGMA_X - L_SIGMAF_LT_full
  diag(SIGMA_U_DIAG) <- 0
  SIGMA_UF_DIAG <- norm(SIGMA_U_DIAG, type = "F")

  return(list(sigma_u = SIGMA_U,
              sigma_u_norm = SIGMA_UF,
              sigma_u_diag = SIGMA_U_DIAG,
              sigma_u_diag_norm = SIGMA_UF_DIAG))
}


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
#' @importFrom utils head
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
