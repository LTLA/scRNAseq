#' Obtain the Segerstolpe pancreas data
#'
#' Download and cache the Segerstolpe pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the pancreas scRNA-seq data from Segerstolpe et al. (2015)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Row data contains fields for the gene symbol and RefSeq transcript IDs corresponding to each gene.
#' Spike-ins are specially labelled with the \code{\link{isSpike}} function.
#'
#' Column metadata were extracted from the \code{Characteristcs} fields of the SDRF file for ArrayExpress E-MTAB-5061.
#' This contains information such as the cell type labels and patient status.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Segerstolpe A et al. (2016). 
#' Single-cell transcriptome profiling of human pancreatic islets in health and type 2 diabetes. 
#' \emph{Cell Metab.} 24(4), 593-607.
#'
#' @examples
#' \dontrun{sce <- SegerstolpePancreasData()}
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
#' @importFrom SummarizedExperiment rowData
SegerstolpePancreasData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("segerstolpe-pancreas", version))
    isSpike(sce, "ERCC") <- grep("^ERCC-[0-9]+", rowData(sce)$refseq)
    sce
}
