#' Greedy algorithm
#' 
#' A greedy version of the greedy peeling algorithm.
#' 
#' The algorithm finds a cluster at each time until the average value in the
#' remainder adjacency matrix is below the mean of the input cluster.
#' 
#' @param Wp Adjacency matrix
#' @param threshold_DSD Threshold value for filtering edges
#' @param lambda Parameter for greedy peeling algorithm
#' @return A list containing the resulting adjacency matrix W_DSD_greedy, 
#'         cluster list Clist, and cluster ID list CID.
#' @export
#' @examples
#' greedy()

greedy <- function(Wp, threshold_DSD, lambda){
  #initialize
  Wp_DSD <- Wp
  Wp_DSD[Wp_DSD < threshold_DSD] <- 0
  
  #store all node lists from each cluster
  nc <- 1
  Clist <- numeric()
  CID <- numeric()
  orig_Node <- 1:length(Wp)
  
  while (length(Clist) < length(Wp) - 2) {
    # do greedy peeling to get the next cluster
    result <- greedy_peeling(Wp_DSD, lambda)
    Clist_temp <- result$Clist
    Node_Seq <- result$Node_Seq
    remaining_node <- result$removing_node
    
    # update adjacency matrix
    Wp_DSD <- Wp_DSD[remaining_node, remaining_node]
    
    # update
    p_rem <- mean(rowMeans(Wp_DSD))
    
    # record nodes
    Clist <- c(Clist, orig_Node[Clist_temp[1:length(Node_Seq)]])
    CID <- c(CID, length(Node_Seq))
    
    # updates
    orig_Node <- orig_Node[remaining_node]
    nc <- nc + 1
  }
  
  # uninformative nodes
  Clist <- c(Clist, orig_Node)
  CID <- c(CID, length(remaining_node))
  
  W_DSD_greedy <- Wp[Clist, Clist]
  
  return(list(W_DSD_greedy = W_DSD_greedy, Clist = Clist, CID = CID))
}
