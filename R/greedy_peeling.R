#' Greedy Peeling
#' 
#' A simple implementation of the greedy algorithm from the "Denser than the Densest subgraph" paper (Tsourakakis 2013)
#' using the generalized objective function (Shuo 2023).
#' 
#' Note that this function only extracts ONE dense subgraph.
#' 
#' @param Wp_DSD Input adjacency matrix
#' @param lambda Parameter for the objective function
#' @return A list with components
#' \itemize{
#'   \item{W_DSD_greedy:}{Reordered adjacency matrix}
#'   \item{Clist:}{Cluster list}
#'   \item{Node_Seq:}{A vector of reordered nodes}
#'   \item{removing_node:}{Nodes to be removed}
#' }
#' @export
#' @examples
#' data(sim)
#' greedy_peeling(sim, 1.4)
#' greedy_peeling(sim, 1.6)

greedy_peeling <- function(Wp_DSD, lambda){
  # Number of nodes in the graph
  N <- nrow(Wp_DSD)
  
  Recording_Matrix <- matrix(0, nrow = N, ncol = 2)
  Recording_Clist <- 1:N
  Wp_temp <- Wp_DSD
  
  for (i in N:1) {
    idxlist_temp <- 1:length(Wp_temp)
    idx_min_temp <- ifelse(length(idxlist_temp) > 1, which.min(colSums(Wp_temp)), 1)
    idxlist_temp <- idxlist_temp[-idx_min_temp]

    if(length(Wp_temp) != 1){
      Wp_temp <- Wp_temp[idxlist_temp, idxlist_temp]
    }
    score_temp <- greedy_objective(Wp_temp, lambda)
    Recording_Matrix[(N - i + 1), ] <- c(Recording_Clist[idx_min_temp], score_temp)
    Recording_Clist <- Recording_Clist[-idx_min_temp]
  }
  
  max_idx <- which.max(Recording_Matrix[, 2])
  removing_node <- Recording_Matrix[1:max_idx, 1]
  Node_Seq <- rev(Recording_Matrix[(max_idx + 1):N, 1])
  Clist <- c(Node_Seq, removing_node)
  W_DSD_greedy <- Wp_DSD[Clist, Clist]
  
  peeling_results <- list(W_DSD_greedy = W_DSD_greedy, 
                          Clist = Clist, 
                          Node_Seq = Node_Seq, 
                          removing_node = removing_node)
  return(peeling_results)
}
