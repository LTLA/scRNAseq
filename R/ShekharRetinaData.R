#' Obtain the Shekhar retina data
#'
#' Download and cache the Shekhar retina single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the mouse retina scRNA-seq data from Shekhar et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' Column metadata contains the cluster identities as reported in the paper.
#' Note that some cells will have \code{NA} identities as they are present in the count matrix but not in the metadata file.
#' These are presumably low-quality cells that were discarded prior to clustering.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Shekhar K et al. (2016). 
#' Comprehensive classification of retinal bipolar neurons by single-cell transcriptomics.
#' \emph{Cell} 166(5), 1308-1323.
#'
#' @examples
#' sce <- ShekharRetinaData()
#' 
#' @export
ShekharRetinaData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("shekhar-retina", version), has.rowdata=FALSE)
    sce
}
