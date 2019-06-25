#' Obtain the Campbell brain data
#'
#' Download and cache the Campbell brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Campbell et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Column metadata is provided in the same form as supplied in GSE93374.
#' This contains information such as the diet of the mice, sex and proposed cell type for each cell.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Campbell R et al. (2017).
#' A molecular census of arcuate hypothalamus and median eminence cell types. 
#' \emph{Nat. Neurosci.} 20, 484-496. 
#'
#' @examples
#' sce <- CampbellBrainData()
#' 
#' @export
CampbellBrainData <- function() {
    version <- "2.0.0"
    .create_sce(file.path("campbell-brain", version), has.rowdata=FALSE)
}
