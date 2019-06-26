#' Obtain the Romanov brain data
#'
#' Download and cache the Romanov brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Romanov et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' Column metadata is provided in the same form as supplied in GSE74672.
#' This contains information such as the reporter gene expressed in each cell, the mouse line, dissection type and so on.
#'
#' Rows corresponding to spike-in transcripts are labelled with the \code{\link{isSpike}} function.
#' Note that some of the spike-in rows have \code{NA} observations for some (but not all) cells.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun,
#' based on code by Vladimir Kiselev and Tallulah Andrews.
#'
#' @references
#' Romanov RA et al. (2017).
#' Molecular interrogation of hypothalamic organization reveals distinct dopamine neuronal subtypes. 
#' \emph{Nat. Neurosci.} 20, 176-188.
#'
#' @examples
#' sce <- RomanovBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
RomanovBrainData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("romanov-brain", version), has.rowdata=FALSE)
    isSpike(sce, "ERCC") <- grep("^ERCC-[0-9]+$", rownames(sce))
    sce
}
