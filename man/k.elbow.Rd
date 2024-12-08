% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kelbow.R
\name{k.elbow}
\alias{k.elbow}
\title{Determine the Elbow Point for k in SCFA}
\usage{
k.elbow(data, CID, Clist, k = length(CID), method = "Sample", epsilon = 1e-06)
}
\arguments{
\item{data}{A data frame or matrix where rows represent observations and columns represent variables.}

\item{CID}{A numeric vector indicating the number of variables associated with each factor.}

\item{Clist}{A numeric vector indicating the column indices of the variables reordered by factors.}

\item{k}{Integer. The maximum number of factors to evaluate. Defaults to the length of \code{CID}.}

\item{method}{Character. The method of residual term estimation: "Sample" or "MLE". Defaults to "Sample".}

\item{epsilon}{Numeric. A small regularization term added to diagonal elements for numerical stability. Defaults to 1e-6.}
}
\value{
A numeric vector where each element represents the residual norm (sum of off-diagonal and diagonal norms)
for the corresponding number of factors (from 0 to \code{k}).
}
\description{
This function computes the residual error for different numbers of factors (k) in a semi-confirmatory
factor analysis (SCFA), helping to identify the optimal number of factors using the "elbow" method.
}
