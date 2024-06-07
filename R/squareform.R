#' Convert a Squareform Matrix to a Vector
#'
#' Given a squareform matrix, this function converts it into a vector containing 
#' the lower triangular part of the matrix (excluding the diagonal).
#'
#' @param mat A square-form matrix.
#' @return A vector containing the lower triangular part of the matrix.
#' @export
#' @examples
#' data(sim)
#' squareform(sim)
#' 
#' mat <- matrix(1:6, nrow = 3)
#' squareform(mat)

squareform <- function(mat) {
  indices <- lower.tri(mat)
  result <- mat[indices]
  return(result)
}