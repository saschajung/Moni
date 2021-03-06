% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateBackground.R
\name{preprocess_ArchS4}
\alias{preprocess_ArchS4}
\title{Preprocess ARCH S4 expression data}
\usage{
preprocess_ArchS4(scriptDir, cor_thr)
}
\arguments{
\item{scriptDir}{Path to a directory containing R-scripts for downloading gene expression data from
ARCH S4.}

\item{cor_thr}{A numerical value above which two samples are considered to be correlated.}
}
\description{
Downloads gene expression data from ArchS4 and selects an uncorrelated subset of samples.
The selected samples are then written into a file.
}
