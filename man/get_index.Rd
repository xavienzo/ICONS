% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{get_index}
\alias{get_index}
\title{Get Indices of Submatrices from a Block-Structured Adjacency Matrix}
\usage{
get_index(blockid, CID, Clist)
}
\arguments{
\item{blockid}{A numeric value or vector indicating the block(s) (submatrix/submatrices) to be selected.
Block IDs should be integers corresponding to positions in the \code{CID} vector.}

\item{CID}{A numeric vector where each element represents the size (number of elements) of each submatrix.}

\item{Clist}{A numeric vector indicating the rearranged index of the full adjacency matrix based on the
submatrices. The length of \code{Clist} should be equal to the sum of the values in \code{CID}.}
}
\value{
A list with two elements:
\describe{
\item{indices}{A numeric vector of the indices in \code{Clist} corresponding to the specified block(s).}
\item{block_ids}{A numeric vector containing the block IDs for each corresponding index in \code{indices}.}
}
}
\description{
This function extracts the indices of the elements corresponding to one or more submatrices
(blocks) in a larger adjacency matrix. The submatrices are defined by a vector containing the
number of elements in each block, and a vector that holds the rearranged index of the full matrix
based on the submatrix structure.
}
\examples{
CID <- c(3, 4, 2)  # 3 nodes in block 1, 4 nodes in block 2, 2 nodes in block 3
Clist <- c(3, 1, 2, 6, 4, 5, 7, 9, 8)  # Rearranged index of the full matrix

# Get indices for block 1 and 2
result <- get_index(blockid = c(1, 2), CID, Clist)
print(result$indices)    # Indices of the rearranged matrix
print(result$block_ids)  # Corresponding block ids
}
