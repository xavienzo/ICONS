#' Greedy Peeling
#' 
#' A simple implementation of the greedy algorithm from the "Denser than the Densest subgraph" paper (Tsourakakis 2013)
#' using the generalized objective function (Shuo 2023).
#' 
#' Note that this function only extracts ONE dense subgraph.
#' 
#' @param Wp_DSD Input adjacency matrix
#' @param lambda Parameter for the objective function
#' @return A list containing:\tabular{ll}{
#'    \code{W_DSD_greedy} \tab Reordered adjacency matrix \cr
#'    \tab \cr
#'    \code{Clist} \tab Cluster list \cr
#'    \tab \cr
#'    \code{Node_Seq} \tab A vector of reordered nodes \cr
#'    \tab \cr
#'    \code{removing_node} \tab Nodes to be removed \cr
#' }
#' @export
#' @examples
#' data(sim)
#' greedy_peeling(sim, 1.4)
#' greedy_peeling(sim, 1.6)

greedy_peeling <- function(Wp_DSD, lambda){
  # Number of nodes in the graph
  N <- nrow(Wp_DSD)
  # N <- sqrt(length(Wp_DSD))
  Recording_Matrix <- matrix(0, nrow = N, ncol = 2)
  Recording_Clist <- 1:N
  Wp_temp <- Wp_DSD
  # if(N == 1){
  #   C <- mean(Wp_temp)
  # }else{
  #   C <- colSums(Wp_temp)
  # }
  C <- colSums(Wp_temp)
  remove_idx <- integer(0)
  ite <- 0
  
  for (i in seq(N, 1, -1)) {
    # Calculate sum ignoring removed indices
    C[remove_idx] <- Inf
    idx_min_temp <- which.min(C)
    # Update remove_idx
    remove_idx <- c(remove_idx, idx_min_temp)
    # Update sum
    C <- C - Wp_temp[idx_min_temp, ]
    C[remove_idx] <- 0
    sum_wp_temp <- sum(C)
    
    # Calculate score_temp
    ite <- ite + 1
    score_temp <- (sum_wp_temp / 2) / ((N - ite) * (2 * lambda))
    Recording_Matrix[N - i + 1, ] <- c(Recording_Clist[idx_min_temp], score_temp)
  }
  
  # Find the index of the maximum score
  max_idx <- which.max(Recording_Matrix[, 2])
  
  # Extract the removing nodes and node sequence
  removing_node <- Recording_Matrix[1:max_idx, 1]
  Node_Seq <- Recording_Matrix[N:(max_idx + 1), 1]
  # Concatenate node sequence and removing nodes
  Clist <- c(Node_Seq, removing_node)
  # Extract subgraph based on Clist
  W_DSD_greedy <- Wp_DSD[Clist, Clist]
  
  return(list(W_DSD_greedy = W_DSD_greedy, 
              Clist = Clist, 
              Node_Seq = Node_Seq, 
              removing_node = removing_node))
}
