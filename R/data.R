#' Simulated data with interconnected community structure
#'
#' A synthetic dataset whose covariance matrix has the interconnected community
#' structure the package targets: several groups of mutually correlated
#' variables, the groups themselves correlated with one another, plus a set of
#' variables that belong to no group.
#'
#' @format A data frame with 100 rows (observations) and 200 numeric columns
#'   (variables).
#'
#' @source Simulated for the ICONS package.
#'
#' @examples
#' data(sim)
#' dim(sim)
#'
#' W <- cor(sim)
#' part <- icons_detect(W, threshold = 0.6)
#' part$sizes
#'
#' plot_matrix(reorder_matrix(W, part), partition = part)
"sim"
