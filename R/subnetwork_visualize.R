#' Visualize a Matrix with ggplot2
#'
#' This function takes a matrix as input and creates a heatmap to visualize it using ggplot2.
#'
#' @param mat A numeric matrix.
#' @return A ggplot2 object representing the heatmap of the matrix.
#' @export
#' @examples
#' mat <- matrix(rnorm(25), nrow = 5)
#' visualize_matrix(mat)
subnetwork_visualize <- function(mat) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }
  
  # Load required packages
  library(ggplot2)
  library(reshape2)
  
  # Convert the matrix to a data frame in long format
  mat_melt <- melt(mat)
  
  # Create the heatmap
  p <- ggplot(mat_melt, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    labs(x = "Column", y = "Row", fill = "Value", title = "Matrix Heatmap") +
    theme_minimal()
  
  return(p)
}