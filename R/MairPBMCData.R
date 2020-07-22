#' Obtain the Mair CITE-seq data
#'
#' Obtain the Mair PBMC targetted CITE-seq data from Mair et al. (2020).
#'
#' @details
#' Column metadata comes directly from the count tables in GSE135325, since there is no separate metadata.
#' This contains information such as donor and batch.
#'
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/mair-pbmc}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Mair C et al. (2020).
#' A Targeted Multi-omic Analysis Approach Measures Protein Expression and 
#' Low-Abundance Transcripts on the Single-Cell Level.
#' \emph{Cell Reports} 31(1), 107499
#'
#' @examples
#' sce <- MairPBMCData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps reducedDim<- 
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment SingleCellExperiment altExp altExpNames
AztekinTailData <- function() {
  version <- "2.4.0"
  
  ##Something like this...
  path <- file.path("mair-pbmc", "2.4.0") #?
  
  rna_counts <- readRDS(file.path(path, "rna_counts.rds"))
  adt_counts <- readRDS(file.path(path, "adt_counts.rds"))
  meta <- readRDS(file.path(path, "coldata.rds"))
  
  sce <- SingleCellExperiment(list(counts = rna_counts))
  colData(sce) <- meta
  altExp(sce) <- SingleCellExperiment(list(counts=adt_counts))
  colData(altExp(sce)) <- meta
  altExpNames(sce) <- "ADT" 
  
  sce
}
