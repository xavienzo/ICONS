#' Simulated data
#'
#' The data contains variables with complex inter-community and intra-community
#' correlations.
#'
#' @format
#' A data frame with 100 subjects and 200 variables.
#' @source null
#' @examples
#' data(sim)
#' matrix <- cor(sim)
#' plotMatrix(data = matrix)
#' results <- dense(matrix, 0.5, 0.5)
#' results$CID

"sim"
