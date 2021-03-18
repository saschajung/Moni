#' Counts to TPM
#'
#' @description Transform a matrix/data frame of gene expression counts into a matrix/data frame of
#' TPM values.
#' @param counts A matrix or data frame object with gene expression counts. Columns correspond to samples and
#' rows to genes.
#' @param len A vector of gene lengths having the same number of entries as there are rows in **counts**,
#' @return A matrix/data frame object with TPM values having the same dimensions as **counts**.
#' @export
tpm3 <- function(counts,len) {
  x <- counts/len
  return(base::t(base::t(x)*1e6/base::colSums(x)))
}

#' Preprocess ARCH S4 expression data
#'
#' @description Downloads gene expression data from ArchS4 and selects an uncorrelated subset of samples.
#' The selected samples are then written into a file.
#' @param scriptDir Path to a directory containing R-scripts for downloading gene expression data from
#' ARCH S4.
#' @param cor_thr A numerical value above which two samples are considered to be correlated.
#' @return NULL
#' @export
preprocess_ArchS4 <- function(scriptDir,cor_thr){
  if(!base::file.exists(base::paste0(scriptDir,"/combinedData.tsv"))){
    files <- base::setdiff(base::list.files(scriptDir,pattern = "*.R"),base::list.files(scriptDir,pattern = "*.tsv"))
    for(f in files){
      #Get the data and run all R scripts first
      base::system(base::paste0("Rscript ","\"",scriptDir,"/",f,"\""))
    }
    files <- base::setdiff(base::list.files(scriptDir,pattern = "*.tsv"),base::c(base::list.files(scriptDir,pattern = "*_ivs.tsv"),"combinedData.tsv"))
    dt_all_test <- data.table::fread(file = base::paste0(scriptDir,"/",files[1]),header = "auto")
    for(f in files[-1]){
      dt <- data.table::fread(file = base::paste0(scriptDir,"/",f),header = "auto")
      dt_all_test <- data.table::merge.data.table(dt_all_test,dt, by = base::c("V1"),all = TRUE)
      #colnames(dt_all) <- paste0("V",1:ncol(dt_all))
    }
    data.table::fwrite(dt_all,file = base::paste0(scriptDir,"/combinedData.tsv"), quote = FALSE, sep= '\t')
  }
  ######## Actual Background Selection #########
  dt_all <- data.table::fread(file = base::paste0(scriptDir,"/combinedData.tsv"),header = FALSE)
  taken <- base::rep(FALSE,base::ncol(dt_all))
  first <- base::sample(base::seq(1,base::ncol(dt_all)),1)
  taken[first] <- TRUE
  dt_set <- dt_all[,base::colnames(dt_all)[first],with = FALSE]
  while(!base::all(taken)){
    colsToTake <- base::setdiff(base::seq(1,base::ncol(dt_all)),base::which(taken))
    s <- base::sample(colsToTake,base::min(base::length(colsToTake),100),replace = FALSE)
    taken[s] <- TRUE
    for(samp in s){
      dt_cor <- base::sapply(1:base::ncol(dt_set),function(x){stats::cor(dt_set[,base::c(x),with = FALSE],dt_all[,base::colnames(dt_all)[samp],with = FALSE])})
      if(base::all(dt_cor < 0.9)){
        dt_set[, (base::colnames(dt_all)[samp]) := dt_all[,base::colnames(dt_all)[samp],with = FALSE]]
      }
    }
  }
  dt_set[, ("V1") := dt_all[,"V1",with = FALSE]]
  data.table::setcolorder(dt_set,base::colnames(dt_set)[base::order(base::colnames(dt_set))])
  data.table::fwrite(dt_set,file = base::paste0(scriptDir,"/Background.tsv"),quote = FALSE, sep = '\t', append = FALSE)
  ######## Actual Background Selection #########
  return(NULL)
}

#' Helper function
#'
#' @description Computes the size of a rectangle above a given empirical cumulative distribution function.
#' @param x A gene expression value.
#' @param empcdf An empirical cumulative distribution function.
#' @param maxVal The maximum expression value to consider. (UNUSED)
#' @return The size of the rectangle.
#' @export
.optLowFun <- function(x,empcdf,maxVal){
  return(x* (1-empcdf(x)))
}

#' Helper function
#'
#' @description Computes the size of a rectangle below a given empirical cumulative distribution function.
#' @param x A gene expression value.
#' @param empcdf An empirical cumulative distribution function.
#' @param maxVal The maximum expression value to consider.
#' @return The size of the rectangle.
#' @export
.optHighFun <- function(x,empcdf,maxVal){
  return(empcdf(x) * (maxVal-x))
}

#' Maximal Rectangle below/above an ecdf
#'
#' @description Computes lower and upper gene expression thresholds by maximizing rectangles below/above a given
#' empirical cumulative distribution function.
#' @param empCDF An empirical cumulative distribution function.
#' @param maxVal The maximum expression value to consider.
#' @return A vector of rectangle sizes.
#' @export
maximizeRectangle <- function(empCDF,maxVal){
  lowThr <- stats::optimize(f = .optLowFun, interval = base::c(0,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum
  highThr <- stats::optimize(f = .optHighFun, interval = base::c(0,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum
  return(base::c(lowThr,highThr))
}

#' Create Threshold Distribution
#'
#' @description Computes a distribution of lower and upper gene expression thresholds of a single gene.
#' @param gexp A vector of gene expression values.
#' @param numBootstrapSamples Number of gene expression samples to draw (with replacement). (DEFAULT: 1000)
#' @return A data frame containing lower and upper gene expression thresholds.
#' @export
createThresholdDist <- function(gexp, numBootstrapSamples = 1000){
  thrs <- base::do.call("rbind",base::lapply(base::seq(1:numBootstrapSamples),function(x){
    samp <- base::sample(gexp,base::length(gexp),replace = TRUE)
    return(Moni::maximizeRectangle(stats::ecdf(samp),base::max(gexp)+1e-4))
  }))
  base::colnames(thrs) <- base::c("Lower","Upper")
  return(thrs)
}

#' Create Threshold Distributions for all genes
#'
#' @description Computes the distribution of lower and upper gene expression thresholds for all genes.
#' @param gexp_mat A matrix of gene expression values with columns corresponding to samples and rows to genes.
#' @param numBootstrapSamples Number of gene expression samples to draw (with replacement). (DEFAULT: 1000)
#' @param norm Normalize each gene by its maximum value? (DEFAULT: FALSE)
#' @return A data frame containing lower and upper gene expression thresholds.
#' @export
createAllThresholdDists <- function(gexp_mat, numBootstrapSamples = 1000, norm = FALSE){
  thr_list <- base::list()
  for(i in 1:base::nrow(gexp_mat)){
    toTest <- base::unname(gexp_mat[base::rownames(gexp_mat)[i],])
    if(norm){
      toTest <- toTest/base::max(toTest)
    }
    thr_list[[base::rownames(gexp_mat)[i]]] <- Moni::createThresholdDist(toTest,numBootstrapSamples)
  }
  return(thr_list)
}

#' Compute P-values for discretizing gene expression values
#'
#' @description Computes p values for a gene to be active or inactive given a threshold distribution.
#' @param expMat A matrix/data frame of gene expression values with columns corresponding to samples and
#' rows to genes.
#' @param thrs_list A list of lower and upper gene expression thresholds per gene.
#' @param p.adj Adjust p values? (DEFAULT: TRUE)
#' @param method Method to use for adjusting p values. Accepted parameter values are the same as for
#' the p.adjust method. (DEFAULT: "fdr")
#' @return A data frame of (adjusted) p values having the same dimension as the input matrix.
#' @export
getPValues <- function(expMat,thrs_list,p.adj = TRUE, method = "fdr"){
  pvals <- base::do.call("rbind",base::lapply(base::rownames(expMat),function(x){
    q <- stats::ecdf(base::c(thrs_list[[x]][,1],thrs_list[[x]][,2]))
    return(q(expMat[x,]))
  }))
  base::rownames(pvals) <- base::rownames(expMat)
  if(p.adj){
    qvals <- base::data.frame(QVal = stats::p.adjust(pvals,method = method))
    base::rownames(qvals) <- base::rownames(expMat)
    return(qvals)
  }else{
    return(pvals)
  }
}

