# Moni - Multiomics data integration unveils core transcriptional regulatory networks governing cell type identity

Reconstruct mechanistic core gene regulatory networks from RNA-seq, epigenetic, protein-protein interaction and TF ChIP-seq datasets.

## Required software
  - bedtools v2.24.2.1 (https://bedtools.readthedocs.io/en/latest/)
  - R v4.0 or greater

## Overview

This R-package provides the basic functionality for network reconstruction. The main functions to be used are **preprocess_ArchS4**, **computeJSDSingle**, **generateNetwork**, **generateNetworkLogic**. If a background gene expression dataset has to be generated, the user first has to download R-scripts from ArchS4 that download the desired files and subsequently call **preprocess_ArchS4** with the directory of the downloaded R-scripts as a parameter. Afterwards, the background can be used to compute the Jennsen-Shannon-Divergence of a query sample with respect to the background for selecting core TFs. Given epigenetic data and a list of core TFs and expressed TFs, **generateNetwork** can then be used to obtaine the interactions between the core and neighboring TFs. Finally, **generateNetworkLogic** will create the logic rules. More documentation about the parameters can be found in the help-files of each function. 

