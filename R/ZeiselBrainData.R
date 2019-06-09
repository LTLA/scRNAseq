#' Obtain the Zeisel brain data
#'
#' Download and cache the Zeisel brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Zeisel et al. (2015)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' Row data contains a single \code{"featureType"} field describing the type of each feature
#' (endogenous genes, mitochondrial genes, spike-in transcripts and repeats).
#' Spike-ins are also specially labelled with the \code{\link{isSpike}} function.
#'
#' Column metadata is provided in the same form as supplied in \url{http://linnarssonlab.org/cortex/}.
#' This contains information such as the cell diameter and the published cell type annotations.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Zeisel A et al. (2015). 
#' Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. 
#' \emph{Science} 347(6226), 1138-42.
#'
#' @examples
#' \dontrun{sce <- ZeiselBrainData()}
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
#' @importFrom SummarizedExperiment rowData
ZeiselBrainData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("zeisel-brain", version))
    isSpike(sce, "ERCC") <- rowData(sce)$featureType=="ERCC"
    sce
}
