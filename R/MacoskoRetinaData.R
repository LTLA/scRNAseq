#' Obtain the Macosko retina data
#'
#' Download and cache the Macosko retina single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the mouse retina scRNA-seq data from Macosko et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' Column metadata contains the cluster identity as reported in the paper, as obtained from the McCarroll website.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Macosko E et al. (2016). 
#' Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. 
#' \emph{Cell} 161(5), 1202-1214.
#'
#' @examples
#' sce <- MacoskoRetinaData()
#' 
#' @export
MacoskoRetinaData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("macosko-retina", version), has.rowdata=FALSE)
    sce
}
