#' Obtain the Chen brain data
#'
#' Download and cache the Chen brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Chen et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Column metadata is provided in the same form as supplied in GSE87544.
#' This contains the putative cell type assigned by the original authors.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Chen R et al. (2017).
#' Single-Cell RNA-Seq reveals hypothalamic cell diversity. 
#' \emph{Cell Rep.} 18, 3227-3241.
#'
#' @examples
#' sce <- ChenBrainData()
#' 
#' @export
ChenBrainData <- function() {
    version <- "2.0.0"
    .create_sce(file.path("chen-brain", version), has.rowdata=FALSE)
}
