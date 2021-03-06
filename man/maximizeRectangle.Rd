% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateBackground.R
\name{maximizeRectangle}
\alias{maximizeRectangle}
\title{Maximal Rectangle below/above an ecdf}
\usage{
maximizeRectangle(empCDF, maxVal)
}
\arguments{
\item{empCDF}{An empirical cumulative distribution function.}

\item{maxVal}{The maximum expression value to consider.}
}
\value{
A vector of rectangle sizes.
}
\description{
Computes lower and upper gene expression thresholds by maximizing rectangles below/above a given
empirical cumulative distribution function.
}
