#' Obtain the Bunis haematopoietic stem and progenitor cell data
#'
#' Obtain the human fetal, newborn, and adult haematopoietic stem and progenitor cell single-cell RNA-seq dataset from Bunis et al. (2021).
#'
#' @param filtered Logical scalar indicating whether to filter out cells that were not used by the authors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata is recreated from GEO using the author-supplied TSV of per-cell annotations, or retrieved from a processed version of the data shared by authors via figshare.
#' This contains information such as the tissue & sample of origin, age group, likely cell type, and Developmental Stage Scoring. Cevelopmental Stage 
#'
#' If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
#' Otherwise, an additional \code{filtered} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors. 
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bacher-tcell}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Daniel Bunis
#'
#' @references
#' Bunis et al. 2021
#'
#' @examples
#' sce <- BunisHSPCData()
#' 
#' @export
BunisHSPCData <- function(filtered=TRUE, rowdata=TRUE) {
    version <- "2.6.0"

    sce <- .create_sce(file.path("bunis-hspc", version))

    if (filtered) {
        sce <- sce[,sce$retained]
        sce$filtered <- NULL
    }
}
