#' Obtain the Xin pancreas data
#'
#' Download and cache the Xin pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the pancreas scRNA-seq data from Xin et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of RPKMs.
#'
#' Row data contains fields for the Entrez ID and symbol for each gene.
#' Column metadata was obtained from the authors (indirectly, via the Hemberg group) 
#' and contains information such as the cell type labels and donor status.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun,
#' using additional metadata obtained by Vladimir Kiselev.
#'
#' @references
#' Xin A et al. (2016).
#' RNA sequencing of single human islet cells reveals type 2 diabetes genes.
#' \emph{Cell Metab.} 24(4), 608-615.
#'
#' @examples
#' sce <- XinPancreasData()
#' 
#' @export
XinPancreasData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("xin-pancreas", version), assays="rpkm")
    sce
}
