#' Objective Function Calculation
#'
#' This function calculates the objective function value based on the input squareform matrix.
#'
#' The objective function is defined as:
#' \deqn{f(C) = \left(\frac{\sum_{i < j} C_{ij}}{n(n-1)/2}\right)^\lambda \cdot \left(\sum_{i < j} C_{ij}\right)^{1-\lambda}}
#'
#' @param C A matrix.
#' @param lambda The lambda parameter used in the objective function (default is 0.5).
#' @return The value of the objective function.
#' @examples
#' C <- matrix(1:6, nrow = 3)
#' greedy_objective(C)
#' @export

greedy_objective <- function(C, lambda) {
  C_vec <- squareform(C)
  value <- ((sum(C_vec)/length(C_vec))^lambda) * ((sum(C_vec))^(1 - lambda))
  return(value)
}