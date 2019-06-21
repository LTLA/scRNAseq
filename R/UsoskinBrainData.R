#' Obtain the Usoskin brain data
#'
#' Download and cache the Usoskin brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Usoskin et al. (2015)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of RPMs.
#'
#' Column metadata is provided in the same form as supplied in External Table 2 of \url{http://linnarssonlab.org/drg/}.
#' This contains information such as the library of origin and the cell type.
#' 
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Usoskin A et al. (2015).
#' Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. 
#' \emph{Nat. Neurosci.} 18(1), 145-53.
#'
#' @examples
#' sce <- UsoskinBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
UsoskinBrainData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("tasic-brain", version), assays="rpms", has.rowdata=FALSE)
    isSpike(sce, "ERCC") <- grep("^ERCC-[0-9]+$", rownames(sce))
    sce
}
