#' Obtain the Bunis haematopoietic stem and progenitor cell data
#'
#' Obtain the human fetal, newborn, and adult haematopoietic stem and progenitor cell single-cell RNA-seq dataset from Bunis et al. (2021).
#'
#' @param filtered Logical scalar or "cells" indicating whether to:
#' \itemize{
#' \item \code{TRUE}: filter out cells that were not used by the authors.
#' \item \code{"cells"}: filter out empty droplets as filtered out by cell ranger.
#' \item \code{FALSE}: no filtering
#' }
#' 
#' @details
#' Column metadata is recreated from GEO using the author-supplied TSV of per-cell annotations, or retrieved from a processed version of the data shared by authors via figshare.
#' This contains information such as the tissue & sample of origin, age group, likely cell type, and Developmental Stage Scoring.
#' Within DevStageScoring element of the column metadata are the applied results ('<cell_type>_scores') of random forest regression trained on the fetal (score = 0) and adult (score = 1) cells of individual cell types indicated by ('<cell_type>_inTraining').
#'
#' If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
#' Otherwise, an additional \code{retained} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bunis-hspc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Daniel Bunis
#'
#' @references
#' Bunis DG et al. (2021). 
#' Single-Cell Mapping of Progressive Fetal-to-Adult Transition in Human Naive T Cells
#' \emph{Cell Rep.} 34(1): 108573
#' 
#' @examples
#' sce <- BunisHSPCData()
#' 
#' @export
BunisHSPCData <- function(filtered=TRUE) {
    version <- "2.6.0"

    sce <- .create_sce(file.path("bunis-hspc", version), has.rowdata = TRUE, has.coldata = FALSE)
    
    hub <- .ExperimentHub()
    colData.path <- file.path("scRNAseq", "bunis-hspc", version, "coldata.rds")
    colData <- hub[hub$rdatapath==colData.path][[1]]

    if (isTRUE(filtered)) {
        keep <- colnames(sce) %in% rownames(colData)[colData$retained]
        sce <- sce[,keep]
        colData$retained <- NULL
    } else if (identical(filtered, "cells")) {
        keep <- colnames(sce) %in% rownames(colData)
        sce <- sce[,keep]
    }

    # Weird performance issue when directly subsetting with rownames.
    # Also, preserve names when filtered=FALSE, though this takes some time.
    m <- match(colnames(sce), rownames(colData))
    colData <- colData[m,, drop = FALSE]
    rownames(colData) <- colnames(sce)
    colData(sce) <- colData

    sce
}
