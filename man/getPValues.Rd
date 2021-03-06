% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateBackground.R
\name{getPValues}
\alias{getPValues}
\title{Compute P-values for discretizing gene expression values}
\usage{
getPValues(expMat, thrs_list, p.adj = TRUE, method = "fdr")
}
\arguments{
\item{expMat}{A matrix/data frame of gene expression values with columns corresponding to samples and
rows to genes.}

\item{thrs_list}{A list of lower and upper gene expression thresholds per gene.}

\item{p.adj}{Adjust p values? (DEFAULT: TRUE)}

\item{method}{Method to use for adjusting p values. Accepted parameter values are the same as for
the p.adjust method. (DEFAULT: "fdr")}
}
\value{
A data frame of (adjusted) p values having the same dimension as the input matrix.
}
\description{
Computes p values for a gene to be active or inactive given a threshold distribution.
}
