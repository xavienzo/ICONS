#' SICERS
#'
#' This function is for parsimonious detector.
#'
#' @param W A n by n matrix of t test -log(p)-values from the raw data, where n denotes the number of nodes.
#' @param r A threshold on the -log(p)-values. We only do the clustering on the significant edges.
#' @param lambda The lambda parameter used in the objective function.
#' @param kmeans_iter Number of iterations for k-means algorithm.
#' @return A list containing:
#'   \describe{
#'     \item{CIDA}{The cluster sizes of every cluster in a power descending order.}
#'     \item{W_SICERS}{Output reordered matrix.}
#'     \item{Clist}{The reordered node index.}
#'   }
#' @examples
#' W <- matrix(1:9, nrow = 3)
#' SICERS_final(W, 0.05, 0.5, 10)
#' @export
SICERS <- function(W, r, lambda, kmeans_iter) {
  # Preprocessing of data
  W1 <- W
  W[W < r] <- 0  # Threshold on the p-values
  z1 <- which(colSums(W) > 0)  # Exclude the isolated nodes
  W <- W[z1, z1]
  degs <- rowSums(W)
  L <- diag(degs) - W  # Laplacian matrix
  
  # Determine the number of clusters K
  U <- eigen(L)$vectors
  U <- U[, ncol(U):1]
  Mk <- matrix(nrow = kmeans_iter, ncol = 2)
  Qual <- matrix(nrow = ncol(U), ncol = kmeans_iter)
  
  for (m in 1:kmeans_iter) {
    Prp_net <- numeric(length = ncol(U))
    for (K in seq(1, ncol(U), by = 5)) {
      C <- kmeans(U[, 1:K], K, iter.max = kmeans_iter)$cluster
      
      indx <- numeric()
      A_net <- numeric()
      net_V <- numeric()
      C_net <- numeric()
      for (k in 1:K) {
        indx <- c(indx, which(C == k))
        net_V[k] <- sum(C == k)
        WC <- W[C == k, C == k]
        C_net[k] <- sum(WC[WC > 0])/2
        A_net[k] <- (net_V[k] * (net_V[k] - 1))/2
      }
      
      Prp_net[K] <- (sum(C_net)^(1 - lambda)) * (sum(C_net)/sum(A_net))^lambda
    }
    
    Kbest <- which(Prp_net == max(Prp_net))
    Kbest <- Kbest[1]  # In case several k's give the same Prp_net value
    Mk[m, ] <- c(Kbest, max(Prp_net))
    
    Qual[, m] <- Prp_net
  }
  
  Kfinal <- Mk[which(Mk[, 2] == max(Mk[, 2])), 1]
  Kfinal <- Kfinal[1]
  
  # Find the cluster ID for each of the nodes
  C <- kmeans(U[, 1:Kfinal], Kfinal, iter.max = kmeans_iter)$cluster
  indx <- numeric()
  A_net <- numeric()
  net_V <- numeric()
  C_net <- numeric()
  for (k in 1:Kfinal) {
    indx <- c(indx, which(C == k))
    net_V[k] <- sum(C == k)
    WC <- W[C == k, C == k]
    C_net[k] <- sum(WC[WC > 0])/2
    A_net[k] <- (net_V[k] * (net_V[k] - 1))/2
  }
  
  diagscore <- (C_net)^(1 - lambda) * (C_net/A_net) * lambda
  diagscore[is.nan(diagscore)] <- 0
  diagscore_sort <- sort(diagscore, decreasing = TRUE)
  diagscore_sortID <- order(diagscore, decreasing = TRUE)
  
  CIDA <- numeric()
  inx_imporance <- numeric()
  for (i in 1:Kfinal) {
    inx_imporance <- c(inx_imporance, which(C == diagscore_sortID[i]))
    CIDA <- c(CIDA, sum(C == diagscore_sortID[i]))
  }
  
  Cindx <- seq_len(nrow(W1))
  Cindx[z1] <- C
  Cindx[setdiff(seq_len(nrow(W1)), z1)] <- -1
  CID <- diagscore_sortID
  Clist <- z1[inx_imporance]
  Clist <- c(Clist, setdiff(seq_len(nrow(W1)), z1))
  CIDA <- c(CIDA, length(setdiff(seq_len(nrow(W1)), z1)))
  W_SICERS <- W1[Clist, Clist]
  
  return(list(CIDA = CIDA, W_SICERS = W_SICERS, Clist = Clist))
}