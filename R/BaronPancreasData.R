#' Obtain the Baron pancreas data
#'
#' Download and cache the Baron pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @param which String specifying the species to get data for.
#'
#' @details
#' This function provides the pancreas scRNA-seq data from Baron et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts for human or mouse.
#'
#' Column metadata is provided in the same form as supplied in GSE84133.
#' This contains information such as the cell type labels and donor ID (for humans) or strain (for mouse).
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Baron M et al. (2017). 
#' Single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. 
#' \emph{Cell Syst.} 3(4), 346-360.
#'
#' @examples
#' sce <- BaronPancreasData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
BaronPancreasData <- function(which=c("human", "mouse")) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("baron-pancreas", version), has.rowdata=FALSE, suffix=which)
    sce
}
