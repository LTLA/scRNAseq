#' Obtain the Tasic brain data
#'
#' Download and cache the Tasic brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Tasic et al. (2015)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Column metadata is provided in the same form as supplied in GSE71585.
#' This contains information such as the reporter gene expressed in each cell, the mouse line, dissection type and so on.
#'
#' Spike-in transcripts are stored as alternative experiments.
#' Note that some of the spike-in rows have \code{NA} observations for some (but not all) cells.
#'
#' The last 9 columns (containing \code{_CTX_} in their names) correspond to no-cell control libraries.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Tasic A et al. (2016). 
#' Adult mouse cortical cell taxonomy revealed by single cell transcriptomics.
#' \emph{Nat. Neurosci.} 19(2), 335-46.
#'
#' @examples
#' sce <- TasicBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
TasicBrainData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("tasic-brain", version), has.rowdata=FALSE)

    status <- ifelse(grepl("^ERCC-[0-9]+$", rownames(sce)), "ERCC", "endogenous")
    splitAltExps(sce, status, ref="endogenous")
}
