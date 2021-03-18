
#' Compute Z-scores
#'
#' @description Computes the Z-score of values given the distribution means and standard deviations.
#' @param x A numeric value or vector of numeric values to standardize.
#' @param m A numeric value or vector of numeric values corresponding to the distribution mean(s).
#' @param s A numeric value or vector of numeric values corresponding to the distribution standard deviation(s).
#' @return A numeric value or vector of numeric values of Z-scores.
#' @export
ztest <- function(x,m,s){
  (x-m)/(s)
}

#' Normalize expression matrix
#'
#' @description Normalizes each gene (row) by its maximum value.
#' @param exp_df A matrix/data frame of gene expression values with rows corresponding to genes
#' and columns to samples.
#' @return A matrix/data frame with normalized expression values.
#' @export
normalizeExpression <- function(exp_df){
  return(base::apply(exp_df,1,function(x){x/base::sum(x)}))
}

#' Compute Regulators of a gene
#'
#' @description Given a data frame of gene regulatory logic rules, computes all regulators.
#' @param network_df A two-column data frame with the first column containing genes and the second
#' column containing the regulatory rules.
#' @return A data frame of unique regulators per gene.
#' @export
getAllRegulators <- function(network_df){
  network_df$Logic <- base::apply(network_df,1,function(x){base::gsub("\\|","",x[2])})
  network_df$Logic <- base::apply(network_df,1,function(x){base::gsub("\\&","",x[2])})
  network_df$Logic <- base::apply(network_df,1,function(x){base::gsub("\\(","",x[2])})
  network_df$Logic <- base::apply(network_df,1,function(x){base::gsub("\\)","",x[2])})
  network_df$Logic <- base::apply(network_df,1,function(x){base::paste(base::unique(base::strsplit(x[2],"\\s+")[[1]]),collapse = ",")})
  network_df$Logic <- base::apply(network_df,1,function(x){if(base::startsWith(x[2],",")){base::substring(x[2],2)}else{x[2]}})

  return(network_df)
}

#' Compute Jennsen-Shannon-Divergence for a single sample
#'
#' @description Compute Jennsen-Shannon-Divergence for a single sample given a data frame of background
#' samples.
#' @param querySample A one-column data frame containing gene expression values. Rownames have to be gene
#' identifiers corresponding to the background data frame.
#' @param background A data frame containing gene expression values. Rownames have to be gene
#' identifiers corresponding to the query sample.
#' @return The value of the Jennsen-Shannon-Divergence.
#' @export
computeJSDSingle <- function(querySample, background){
  merged_df <- base::merge(querySample,background,by="row.names")
  base::rownames(merged_df) <- merged_df$Row.names
  merged_df$Row.names <- NULL
  cors <- base::apply(merged_df,2,function(x){stats::cor(merged_df[,1],x)})
  merged_df <- merged_df[,base::c(1,base::which(cors <= 0.7))]
  merged_df <- base::t(Moni::normalizeExpression(merged_df))

  ideal_distr_vector <- base::rep(0,base::ncol(merged_df))
  ideal_distr_vector[1] = 1

  JSD_based_expr_score <- base::apply(merged_df,1, function(x) {
    return(base::unname(philentropy::JSD(base::rbind(x,ideal_distr_vector))))
  })

  return(JSD_based_expr_score)
}

#' Compute Jennsen-Shannon-Divergence for all samples
#'
#' @description Compute Jennsen-Shannon-Divergence for a all samples in a set of background
#' samples.
#' @param background A data frame containing gene expression values. Rownames have to be gene
#' identifiers corresponding to the query sample.
#' @return A data frame of Jennsen-Shannon-Divergence values and their ranks with respect to all samples.
#' @export
computeJSDBackground <- function(background){
  jsd <- base::matrix(data = -1,nrow = base::nrow(background), ncol = base::ncol(background))
  for(i in 1:base::ncol(background)){
    base::suppressMessages(jsd[,i] <- Moni::computeJSDSingle(background[,i,drop = FALSE],background[,-i]))
  }
  base::rownames(jsd) <- base::rownames(background)
  base::colnames(jsd) <- base::colnames(background)
  jsd_ranks <- base::apply(jsd,2,rank)
  base::rownames(jsd_ranks) <- base::rownames(background)
  base::colnames(jsd_ranks) <- base::colnames(background)
  return(base::list(JSD = jsd, JSD_Ranks = jsd_ranks))
}

#' Compute inclusion list
#'
#' @description Identifies a set of genes whose Jennsen-Shannon-Divergence value could be considered
#' significant.
#' @param query_sample A one-column data frame containing gene expression values. Rownames have to be gene
#' identifiers corresponding to the background data frame.
#' @param jsd_ranks The ranked Jennsen-Shannon-Divergence values for all samples in the background.
#' @param background A data frame containing gene expression values. Rownames have to be gene
#' identifiers corresponding to the query sample.
#' @param z_score_cutoff A Z-score below which a gene is considered significant. (DEFAULT: -1.5)
#' @return A vector of genes considered significantly unique.
#' @export
computeInclusionList <- function(query_sample, jsd_ranks,background,z_score_cutoff = -1.5){
  commongenes <- base::intersect(base::rownames(query_sample),base::intersect(base::rownames(jsd_ranks),base::rownames(background)))
  query_sample <- query_sample[commongenes,,drop = FALSE]
  jsd_ranks <- jsd_ranks[commongenes,]
  background <- background[commongenes,]
  query_sample <- query_sample[base::rownames(background),,drop = FALSE]
  base::suppressMessages(query_jsd <- Moni::computeJSDSingle(query_sample,background))
  query_jsd_ranks <- base::rank(query_jsd)
  cors <- base::apply(background,2,function(x){stats::cor(query_sample,x)})
  query_jsd_ranks <- query_jsd_ranks[base::rownames(jsd_ranks)]
  incl_samples <- base::colnames(jsd_ranks)[base::which(cors <= 0.7)]
  incl_samples <- base::which(base::colnames(jsd_ranks) %in% incl_samples)
  z <- Moni::ztest(query_jsd_ranks,base::sapply(base::rownames(query_sample),function(x){base::mean(jsd_ranks[x,incl_samples])}),base::sapply(base::rownames(query_sample),function(x){stats::sd(jsd_ranks[x,incl_samples])}))
  incl <- base::names(z)[base::which(z <= z_score_cutoff)]
  return(incl)
}

#' Compute logic rules of genes
#'
#' @description Puts together the gene regulatory rules of all genes to be considered in a network
#' @param df A data frame with TF binding sites in regulatory regions of genes.
#' @param ppi A data frame of protein-protein interactions.
#' @param directed Should the network of overlapping binding sites be considered directed or undirected?
#' @param mode Either "weak" or "strong". Should weakly or strongly connected componentes of overlapping binding
#' sites be computed?
#' @return A data frame of gene regulatory logic rules for each gene.
#' @export
writeLogic <- function(df,ppi,directed,mode){
  targetGenes <- base::unique(df$V9)
  logicRegulation_df <- base::data.frame(Gene = targetGenes, Logic = base::rep("",base::length(targetGenes)), stringsAsFactors = FALSE)
  for(i in 1:base::length(targetGenes))
  {
    logic <- base::c()
    ints <- df[base::which(df$V9 == targetGenes[i]),base::c(4,16)]
    graph <- igraph::graph_from_data_frame(ints,directed = directed)
    SCC <- igraph::clusters(graph,mode=mode)
    clusts <- SCC$membership
    base::names(clusts) <- base::gsub("_[0-9]+","", base::names(clusts))
    clusts <- base::cbind.data.frame(base::names(clusts),clusts,stringsAsFactors=FALSE)
    clusts <- clusts[!base::duplicated(clusts),]
    for(c in 1:base::length(SCC$csize))
    {
      cur_clust <- clusts[clusts[,2] == c,1]
      if(base::length(cur_clust) == 1){
        logic <- base::c(logic,cur_clust[1])
      }else{
        ppi_subs <- ppi[ppi$symbol1 %in% cur_clust & ppi$symbol2 %in% cur_clust,]
        ppi_graph <- igraph::graph_from_data_frame(ppi_subs[,1:2],directed = FALSE)
        ppi_graph_comp <- igraph::clusters(ppi_graph)
        if(base::length(ppi_graph_comp$csize) > 0){
          for(ppi_c in 1:base::length(ppi_graph_comp$csize)){
            if(ppi_graph_comp$csize[ppi_c] == 1){
              logic <- base::c(logic,base::names(ppi_graph_comp$membership)[base::which(ppi_graph_comp$membership == ppi_c)])
            }else{
              logic <- base::c(logic,base::paste0("( ",base::paste(base::names(ppi_graph_comp$membership)[base::which(ppi_graph_comp$membership == ppi_c)], collapse = " & "), " )"))
            }
          }
        }

        remaining <- base::setdiff(cur_clust,base::names(ppi_graph_comp$membership))
        if(base::length(remaining) > 0){
          for(j in 1:base::length(remaining)){
            logic <- base::c(logic, remaining[j])
          }
        }
      }
    }
    regMotifs <- df[base::which(df$V9 == targetGenes[i]),4]
    regMotifs <- base::setdiff(regMotifs,base::names(SCC$membership))
    regMotifs <- base::gsub("_[0-9]+","", regMotifs)
    regMotifs <- regMotifs[!base::duplicated(regMotifs)]
    if(base::length(regMotifs) > 0){
      logic <- base::c(logic, regMotifs)
    }

    logicRegulation_df[base::which(logicRegulation_df$Gene == targetGenes[i]),2] <- base::paste(base::unique(logic),collapse = " | ")

  }
  return(logicRegulation_df)
}

#' Generate Network Logics
#'
#' @description Puts together the gene regulatory rules of all genes to be considered in a network
#' @param data_prom A data frame with TF binding sites in (active) promoter regions.
#' @param data_enh A data frame with TF binding sites in (active) enhancer regions.
#' @param ppi A data frame of protein-protein interactions.
#' @param isStrict Should weakly or strongly connected componentes of overlapping binding
#' sites be computed?
#' @return A list containing only the promoter/enhancer network rules, the combined network rules
#' and the unique regulators per gene in enhancers and promoters.
#' @export
generateNetworkLogic <- function(data_prom,data_enh,ppi, isStrict = F){

  base::colnames(ppi) <- base::c("symbol1","symbol2","Score")

  idx <- base::which(!(data_prom$V9 %in% base::union(base::unique(data_prom$V4),base::unique(data_enh$V4))))
  while(base::length(idx) > 0){
    data_prom<- data_prom[-idx,]
    idx <- base::which(!(data_prom$V9 %in% base::union(base::unique(data_prom$V4),base::unique(data_enh$V4))))
  }
  idx <- base::which(!(data_enh$V12 %in% base::union(base::unique(data_prom$V4),base::unique(data_enh$V4))))
  while(base::length(idx) > 0){
    data_enh <- data_enh[-idx,]
    idx <- base::which(!(data_enh$V12 %in% base::union(base::unique(data_prom$V4),base::unique(data_enh$V4))))
  }

  if(base::nrow(data_prom) == 0 | base::nrow(data_enh) == 0){
    return(base::list())
  }

  data_prom <- data_prom[base::order(data_prom$V6,data_prom$V7,data_prom$V8),]
  data_enh <- data_enh[base::order(data_enh$V6,data_enh$V7,data_enh$V8),]

  uniqueRegulators <- data_prom[,1:4]
  uniqueRegulators <- uniqueRegulators[!base::duplicated(uniqueRegulators),]
  uniqueRegulators <- base::cbind.data.frame(uniqueRegulators,base::seq(1,base::nrow(uniqueRegulators)))
  base::colnames(uniqueRegulators)[5] <- "V5"
  uniqueRegulators$V5 <- base::paste0(uniqueRegulators$V4,"_",uniqueRegulators$V5)

  data_prom <-plyr::join(data_prom,uniqueRegulators, by = base::c("V1","V2","V3","V4"))
  data_prom$V4 <- data_prom[,13]
  data_prom <- data_prom[,-13]

  relations_prom <- bedr::bedr(engine = "bedtools",
                          params = base::c("-loj -r -f 0.62"),
                          input = base::list(a = data_prom, b = data_prom),
                          method = "intersect",
                          tmpDir = "./",
                          deleteTmpDir = TRUE,
                          verbose = TRUE)
  relations_prom$V25 <- base::abs(relations_prom$V5-relations_prom$V17)
  targetGenes <- base::unique(relations_prom$V9)

  uniqueRegulators_enh <- data_enh[,1:4]
  uniqueRegulators_enh <- uniqueRegulators_enh[!base::duplicated(uniqueRegulators_enh),]
  uniqueRegulators_enh <- base::cbind.data.frame(uniqueRegulators_enh,base::seq(1,base::nrow(uniqueRegulators_enh)))
  base::colnames(uniqueRegulators_enh)[5] <- "V5"
  uniqueRegulators_enh$V5 <- base::paste0(uniqueRegulators_enh$V4,"_",uniqueRegulators_enh$V5)

  data_enh <-plyr::join(data_enh,uniqueRegulators_enh, by = base::c("V1","V2","V3","V4"))
  data_enh$V4 <- data_enh[,15]
  data_enh <- data_enh[,-15]

  relations_enh <- bedr::bedr(engine = "bedtools",
             params = base::c("-loj -r -f 0.62"),
             input = base::list(a = data_enh, b = data_enh),
             method = "intersect",
             tmpDir = "./",
             deleteTmpDir = TRUE,
             verbose = TRUE)
  targetGenes <- base::unique(relations_enh$V12)
  relations_enh <- relations_enh[,-base::c(9,10,11,23,25)]
  relations_enh <- base::cbind.data.frame(relations_enh[,1:9],base::data.frame(A = base::rep(1,base::nrow(relations_enh))),relations_enh[,10:21],base::data.frame(B = base::rep(1,base::nrow(relations_enh))),relations_enh[,22:23],stringsAsFactors = FALSE)
  relations_enh$V25 <- base::abs(relations_enh$V5-relations_enh$V19)
  base::colnames(relations_enh) <- base::colnames(relations_prom)

  logicRegulation_df <- Moni::writeLogic(relations_prom,ppi,directed = TRUE, mode = base::ifelse(isStrict,"strong","weak"))
  logicRegulation_enh_df <- Moni::writeLogic(relations_enh,ppi,directed = TRUE, mode = base::ifelse(isStrict,"strong","weak"))
  base::colnames(logicRegulation_enh_df) <- base::c("Gene","EnhLogic")

  logics <- plyr::join(logicRegulation_df,logicRegulation_enh_df,by = base::c("Gene"), type = "left")
  logics$EnhLogic[base::which(base::is.na(logics$EnhLogic))] <- "FALSE"
  logics$Logic[base::which(base::is.na(logics$Logic))] <- "FALSE"
  promNetwork <- logics[,base::c(1,2)]
  enhNetwork <- logics[,base::c(1,3)]
  logics$Logic <- base::apply(logics,1,function(x){base::paste0("( ",x[2]," ) & ( ",x[3], " )")})
  logics <- logics[,-3]
  targetGenes <- logics$Gene
  base::colnames(logics) <- base::c("Gene","Logic")

  selfRegs <- targetGenes[base::sapply(targetGenes,function(x){
    tmp <- logics$Gene[base::grepl(x,logics$Logic)]
    return(base::all(tmp == x))
  })]
  if(base::length(selfRegs) > 0){
    logics <- logics[logics$Gene %in% base::setdiff(targetGenes,selfRegs),]
  }
  if(base::nrow(logics) == 0){
    return(base::list())
  }
  promNetwork <- promNetwork[promNetwork[,1] %in% logics$Gene, ]
  enhNetwork <- enhNetwork[enhNetwork[,1] %in% logics$Gene, ]

  promRegulators <-Moni::getAllRegulators(promNetwork)
  enhRegulators <-Moni::getAllRegulators(enhNetwork)

  return(base::list(Logic = logics, PromoterNetwork = promNetwork, EnhancerNetwork = enhNetwork, PromoterRegulators = promRegulators, EnhancerRegulators = enhRegulators))
}

#' Generate networks from .bed files
#'
#' @description Puts together the gene regulatory rules of all genes to be considered in a network
#' @param prom_df A data frame in .bed format. Each line must be a H3K4me3 peak.
#' @param enh_df A data frame in .bed format. Each line must be a H3K27ac peak.
#' @param exprTFs_df A one column data frame of all expressed genes.
#' @param Dnase_df A data frame in .bed format. Each line must be a DNase-seq peak
#' @param tfs_df A one column data frame of core TFs to consider.
#' @param Cistrome_ChIPseq_sorted_df A data frame of TF ChIP-seq binding sites. Columns must correspond to:
#' Chromosome, Start, End, TF (in this order).
#' @param Processed_GeneHancer_df Data frame of the processed gene hancer data.
#' @param KnownPromoters_df Data frame of known promoter regions (from AnimalTFDB).
#' @return A list of data frames corresponding to TF binding sites in active promoters and enhancers for core
#' TFs and non-core mutually connected neighbors.
#' @export
generateNetwork <- function(prom_df, enh_df, exprTFs_df, Dnase_df, tfs_df, Cistrome_ChIPseq_sorted_df, Processed_GeneHancer_df, KnownPromoters_df){

  exprTFs_df <- exprTFs_df[base::order(exprTFs_df[,1]),,drop = FALSE]
  coreTFs_df <- tfs_df[base::order(tfs_df[,1]),,drop = FALSE]
  nonCoreTFs_df <- exprTFs_df[base::which(!(exprTFs_df[,1] %in% coreTFs_sorted[,1])),,drop=FALSE]

  Cistrome_ChIPseq_exprTFs_df <- Cistrome_ChIPseq_sorted_df[base::which(Cistrome_ChIPseq_sorted_df[,4] %in% exprTFs_df[,1]),1:5,drop = FALSE]
  Cistrome_ChIPseq_coreTFs_df <- Cistrome_ChIPseq_sorted_df[base::which(Cistrome_ChIPseq_exprTFs_df[,4] %in% coreTFs_df[,1]),1:5,drop = FALSE]
  Cistrome_ChIPseq_nonCoreTFs_df <- Cistrome_ChIPseq_sorted_df[base::which(Cistrome_ChIPseq_exprTFs_df[,4] %in% nonCoreTFs_df[,1]),1:5,drop = FALSE]

  Processed_GeneHancer_df <- Processed_GeneHancer_df[base::order(Processed_GeneHancer_df[,7]),,drop = FALSE]
  Processed_GeneHancer_df <- Processed_GeneHancer_df[base::which(Processed_GeneHancer_df[,8] > 10),,drop = FALSE]

  potentialEnhancers_df <- Processed_GeneHancer_df[base::which(Processed_GeneHancer_df[,7] %in% coreTFs_df[,1]),1:8,drop = FALSE]

  activeEnhancers_df <- bedr::bedr(engine = "bedtools",
                                  params = base::c(""),
                                  input = base::list(a = potentialEnhancers_df, b = enh_df),
                                  method = "intersect",
                                  tmpDir = "./",
                                  deleteTmpDir = TRUE,
                                  verbose = TRUE
                        )

  coreTFs_in_ActiveEnhancers_df <- bedr::bedr(engine = "bedtools",
                                              params = base::c("-wo -F 1"),
                                              input = base::list(a = activeEnhancers_df, b = Cistrome_ChIPseq_coreTFs_df),
                                              method = "intersect",
                                              tmpDir = "./",
                                              deleteTmpDir = TRUE,
                                              verbose = TRUE
                                    )

  nonCoreTFs_in_ActiveEnhancers_df <- bedr::bedr(engine = "bedtools",
                                              params = base::c("-wo -F 1"),
                                              input = base::list(a = activeEnhancers_df, b = Cistrome_ChIPseq_nonCoreTFs_df),
                                              method = "intersect",
                                              tmpDir = "./",
                                              deleteTmpDir = TRUE,
                                              verbose = TRUE
                                  )

  tfs_in_ActiveEnhancers_df <- bedr::bedr(engine = "bedtools",
                                                 params = base::c("-wo -F 1"),
                                                 input = base::list(a = activeEnhancers_df, b = Cistrome_ChIPseq_exprTFs_df),
                                                 method = "intersect",
                                                 tmpDir = "./",
                                                 deleteTmpDir = TRUE,
                                                 verbose = TRUE
                               )
  tfs_in_ActiveEnhancers_df <- tfs_in_ActiveEnhancers_df[base::order(tfs_in_ActiveEnhancers_df[,12]),,drop = FALSE]

  potentialEnhancers_neigh_df <- Processed_GeneHancer_df[base::which(Processed_GeneHancer_df[,7] %in% tfs_in_ActiveEnhancers_df[,12]),,drop = FALSE]

  activeEnhancers_neigh_df <- bedr::bedr(engine = "bedtools",
                                   params = base::c(""),
                                   input = base::list(a = potentialEnhancers_neigh_df, b = enh_df),
                                   method = "intersect",
                                   tmpDir = "./",
                                   deleteTmpDir = TRUE,
                                   verbose = TRUE
  )

  coreTFs_in_ActiveEnhancers_neigh_df <- bedr::bedr(engine = "bedtools",
                                              params = base::c("-wo -F 1"),
                                              input = base::list(a = activeEnhancers_neigh_df, b = Cistrome_ChIPseq_coreTFs_df),
                                              method = "intersect",
                                              tmpDir = "./",
                                              deleteTmpDir = TRUE,
                                              verbose = TRUE
  )

  tfs_in_ActiveEnhancers_neigh_df <- bedr::bedr(engine = "bedtools",
                                          params = base::c("-wo -F 1"),
                                          input = base::list(a = activeEnhancers_neigh_df, b = Cistrome_ChIPseq_exprTFs_df),
                                          method = "intersect",
                                          tmpDir = "./",
                                          deleteTmpDir = TRUE,
                                          verbose = TRUE
  )
  tfs_in_ActiveEnhancers_neigh_df <- tfs_in_ActiveEnhancers_neigh_df[base::order(tfs_in_ActiveEnhancers_neigh_df[,12]),,drop = FALSE]

  reg_enh_df <- base::rbind.data.frame(coreTFs_df,coreTFs_in_ActiveEnhancers_neigh_df[,7,drop=FALSE])
  reg_enh_df <- reg_enh_df[!base::duplicated(reg_enh_df),,drop=FALSE]
  reg_enh_df <- reg_enh_df[base::order(reg_enh_df[,1]),,drop=FALSE]

  tmp <- base::rbind.data.frame(tfs_in_ActiveEnhancers_df,tfs_in_ActiveEnhancers_neigh_df)
  tmp <- tmp[base::order(tmp[,12]),]

  Enhancers <- tmp[base::which(tmp[,12] %in% reg_enh_df[,1]),base::c(9:13,1:8,14)]
  Enhancers <- Enhancers[!base::duplicated(Enhancers),]
  Enhancers <- Enhancers[base::order(Enhancers[,1],Enhancers[,2],Enhancers[,3]),]

  Enhancers_Acc <- bedr::bedr(engine = "bedtools",
                              params = base::c("-e -wa -f 0.50 -F 0.50"),
                              input = base::list(a = Enhancers, b = Dnase_df),
                              method = "intersect",
                              tmpDir = "./",
                              deleteTmpDir = TRUE,
                              verbose = TRUE
                    )


  # Promoters
  KnownPromoters_df <- KnownPromoters_df[base::order(KnownPromoters_df[,4]),,drop=FALSE]

  potentialPromoters_df <- KnownPromoters_df[base::which(KnownPromoters_df[,4] %in% coreTFs_df[,1]),1:6,drop=FALSE]

  activePromoters_df <- bedr::bedr(engine = "bedtools",
                                   params = base::c(""),
                                   input = base::list(a = potentialPromoters_df, b = prom_df),
                                   method = "intersect",
                                   tmpDir = "./",
                                   deleteTmpDir = TRUE,
                                   verbose = TRUE
  )

  coreTFs_in_ActivePromoters_df <- bedr::bedr(engine = "bedtools",
                                              params = base::c("-wo -F 1"),
                                              input = base::list(a = activePromoters_df, b = Cistrome_ChIPseq_coreTFs_df),
                                              method = "intersect",
                                              tmpDir = "./",
                                              deleteTmpDir = TRUE,
                                              verbose = TRUE
  )

  nonCoreTFs_in_ActivePromoters_df <- bedr::bedr(engine = "bedtools",
                                                 params = base::c("-wo -F 1"),
                                                 input = base::list(a = activePromoters_df, b = Cistrome_ChIPseq_nonCoreTFs_df),
                                                 method = "intersect",
                                                 tmpDir = "./",
                                                 deleteTmpDir = TRUE,
                                                 verbose = TRUE
  )

  tfs_in_ActivePromoters_df <- bedr::bedr(engine = "bedtools",
                                          params = base::c("-wo -F 1"),
                                          input = base::list(a = activePromoters_df, b = Cistrome_ChIPseq_exprTFs_df),
                                          method = "intersect",
                                          tmpDir = "./",
                                          deleteTmpDir = TRUE,
                                          verbose = TRUE
  )

  tfs_in_ActivePromoters_df <- tfs_in_ActivePromoters_df[base::order(tfs_in_ActivePromoters_df[,10]),]

  potentialPromoters_neigh_df <- KnownPromoters_df[base::which(KnownPromoters_df[,4] %in% tfs_in_ActivePromoters_df[,10]),1:6]

  activePromoters_neigh_df <- bedr::bedr(engine = "bedtools",
                                         params = base::c(""),
                                         input = base::list(a = potentialPromoters_neigh_df, b = prom_df),
                                         method = "intersect",
                                         tmpDir = "./",
                                         deleteTmpDir = TRUE,
                                         verbose = TRUE
  )

  coreTFs_in_ActivePromoters_neigh_df <- bedr::bedr(engine = "bedtools",
                                                    params = base::c("-wo -F 1"),
                                                    input = base::list(a = activePromoters_neigh_df, b = Cistrome_ChIPseq_coreTFs_df),
                                                    method = "intersect",
                                                    tmpDir = "./",
                                                    deleteTmpDir = TRUE,
                                                    verbose = TRUE
  )

  tfs_in_ActivePromoters_neigh_df <- bedr::bedr(engine = "bedtools",
                                                params = base::c("-wo -F 1"),
                                                input = base::list(a = activePromoters_neigh_df, b = Cistrome_ChIPseq_exprTFs_df),
                                                method = "intersect",
                                                tmpDir = "./",
                                                deleteTmpDir = TRUE,
                                                verbose = TRUE
  )

  reg_prom_df <- base::rbind.data.frame(coreTFs_df,coreTFs_in_ActivePromoters_neigh_df[,4,drop=FALSE])
  reg_prom_df <- reg_prom_df[!base::duplicated(reg_prom_df),,drop=FALSE]
  reg_prom_df <- reg_prom_df[base::order(reg_prom_df[,1]),,drop=FALSE]

  tmp <- base::rbind.data.frame(tfs_in_ActivePromoters_df,tfs_in_ActivePromoters_neigh_df)
  tmp <- tmp[base::order(tmp[,10]),]

  Promoters <- tmp[base::which(tmp[,10] %in% reg_prom_df[,1]),base::c(7:11,1:6,12)]
  Promoters <- Promoters[!base::duplicated(Promoters),]
  Promoters <- Promoters[base::order(Promoters[,1],Promoters[,2],Promoters[,3]),]

  Promoters_Acc <- bedr::bedr(engine = "bedtools",
                              params = base::c("-e -wa -f 0.50 -F 0.50"),
                              input = base::list(a = Promoters, b = Dnase_df),
                              method = "intersect",
                              tmpDir = "./",
                              deleteTmpDir = TRUE,
                              verbose = TRUE
  )

  return(base::list(Promoters = Promoters_Acc, Enhancers = Enhancers_Acc))
}



