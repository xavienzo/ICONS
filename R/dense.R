#' Greedy peeling algorithm with \eqn{l_0}{l_0} shrinkage
#'
#' A simple implementation of the greedy algorithm from the "Denser than the Densest subgraph" paper (Tsourakakis, 2013)
#' using the generalized objective function (Chen, 2023).
#'
#' Note that this function only extracts ONE dense subgraph.
#'
#' @param W_greedy Input adjacency matrix
#' @param lambda Tuning parameter \eqn{\lambda}{lambda} for greedy peeling algorithm
#' @return A list containing:\tabular{ll}{
#'    \code{W_interim} \tab Reordered adjacency matrix \cr
#'    \tab \cr
#'    \code{Clist} \tab Node index of the reordered adjacency matrix \cr
#'    \tab \cr
#'    \code{Node_Seq} \tab Node index of the extracted dense subgraph \cr
#'    \tab \cr
#'    \code{removing_node} \tab Index of nodes not in the extracted dense subgraph \cr
#' }
#' @references
#' Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
#' Identifying covariate-related subnetworks for whole-brain connectome analysis. Biostatistics (Oxford, England), 25(2), 541–558. https://doi.org/10.1093/biostatistics/kxad007
#'
#' Tsourakakis, C., Bonchi, F., Gionis, A., Gullo, F., & Tsiarli, M. (2013). Denser than the densest subgraph: Extracting optimal quasi-cliques with quality guarantees.
#' Proceedings of the 19th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining.
#' @export
#' @examples
#' data(sim)
#' matrix <- cor(sim)
#' greedy_peeling(matrix, 0.6)
#' greedy_peeling(matrix, 0.3)

greedy_peeling <- function(W_greedy, lambda){
  # Number of nodes in the graph
  N <- nrow(W_greedy)
  # N <- sqrt(length(W_greedy))
  Recording_Matrix <- matrix(0, nrow = N, ncol = 2)
  Recording_Clist <- 1:N
  W_temp <- W_greedy

  C <- colSums(W_temp)
  remove_idx <- integer(0)
  ite <- 0

  for (i in seq(N, 1, -1)) {
    # Calculate sum ignoring removed indices
    C[remove_idx] <- Inf
    idx_min_temp <- which.min(C)
    # Update remove_idx
    remove_idx <- c(remove_idx, idx_min_temp)
    # Update sum
    C <- C - W_temp[idx_min_temp, ]
    C[remove_idx] <- 0
    sum_W_temp <- sum(C)

    # Calculate score_temp
    ite <- ite + 1
    score_temp <- (sum_W_temp / 2) / ((N - ite) ^ (2 * lambda))
    Recording_Matrix[N - i + 1, ] <- c(Recording_Clist[idx_min_temp], score_temp)
  }

  # Find the index of the maximum score
  max_idx <- which.max(Recording_Matrix[, 2])

  # Extract the removing nodes and node sequence
  removing_node <- Recording_Matrix[1:max_idx, 1]
  if (length(removing_node) != 1) {
    Node_Seq <- Recording_Matrix[N:(max_idx + 1), 1]
  } else {
    Node_Seq <- Recording_Matrix[, 1]
  }
  # Concatenate node sequence and removing nodes
  Clist <- c(Node_Seq, removing_node)
  # Extract subgraph based on Clist
  W_interim <- W_greedy[Clist, Clist]

  return(list(W_interim = W_interim,
              Clist = Clist,
              Node_Seq = Node_Seq,
              removing_node = removing_node))
}


#' Adaptive dense subgraph extraction
#'
#' Extract dense subgraphs from the entire graph using the greedy algorithm with \eqn{l_0}{l_0} shrinkage.
#'
#' @param W_original Input adjacency matrix
#' @param threshold Threshold value for filtering edges, default to 0.5
#' @param lambda Tuning parameter \eqn{\lambda}{lambda} for greedy peeling algorithm, default to 0.5
#' @return A list containing:\tabular{ll}{
#'    \code{W_dense} \tab Reordered adjacency matrix \cr
#'    \tab \cr
#'    \code{Clist} \tab Node index of the reordered adjacency matrix with detected dense subgraphs\cr
#'    \tab \cr
#'    \code{CID} \tab Number of nodes in each dense subgraph \cr
#' }
#' @references
#' Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
#' Identifying covariate-related subnetworks for whole-brain connectome analysis.
#' Biostatistics (Oxford, England), 25(2), 541–558. https://doi.org/10.1093/biostatistics/kxad007
#'
#' Wu, Q., Huang, X., Culbreth, A. J., Waltz, J. A., Hong, L. E., & Chen, S. (2022). Extracting brain disease‐related connectome subgraphs by adaptive dense subgraph discovery.
#' Biometrics, 78(4), 1566–1578. https://doi.org/10.1111/biom.13537
#' @export
#' @examples
#' data(sim)
#' matrix <- cor(sim)
#' dense(matrix, 0.6, 0.5)

dense <- function(W_original, threshold = 0.5, lambda = 0.5){

  if (!all(diag(W_original) == 0)) {
    W_original <- W_original - diag(diag(W_original))
  }
  #initialize
  W_greedy <- W_original
  W_greedy[W_greedy < threshold] <- 0

  #store all node lists from each cluster
  Clist <- numeric()
  CID <- numeric()
  orig_Node <- 1:length(W_original)

  while (length(Clist) < ncol(W_original) - 1) {
    # perform greedy peeling to detect a dense subgraph, remove it, perform on the remaining nodes
    result <- greedy_peeling(W_greedy, lambda)
    Clist_temp <- result$Clist
    Node_Seq <- result$Node_Seq
    remaining_node <- result$removing_node

    # record nodes
    Clist <- c(Clist, orig_Node[Clist_temp[1:length(Node_Seq)]])
    CID <- c(CID, length(Node_Seq))

    if (length(remaining_node) == 1) break
    # updates
    orig_Node <- orig_Node[remaining_node]
    W_greedy <- W_greedy[remaining_node, remaining_node]
  }

  W_dense <- W_original[Clist, Clist]

  return(list(W_dense = W_dense, Clist = Clist, CID = CID))
}


