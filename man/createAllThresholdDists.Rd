% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateBackground.R
\name{createAllThresholdDists}
\alias{createAllThresholdDists}
\title{Create Threshold Distributions for all genes}
\usage{
createAllThresholdDists(gexp_mat, numBootstrapSamples = 1000, norm = FALSE)
}
\arguments{
\item{gexp_mat}{A matrix of gene expression values with columns corresponding to samples and rows to genes.}

\item{numBootstrapSamples}{Number of gene expression samples to draw (with replacement). (DEFAULT: 1000)}

\item{norm}{Normalize each gene by its maximum value? (DEFAULT: FALSE)}
}
\value{
A data frame containing lower and upper gene expression thresholds.
}
\description{
Computes the distribution of lower and upper gene expression thresholds for all genes.
}
