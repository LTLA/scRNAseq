#' Obtain the Pollen radial glia data
#'
#' Obtain the human radial glia single-cell RNA-seq dataset from Pollen et al. (2017).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata includes the anatomical source, sample of origin, presumed cell type and assorted alignment statistics.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' This is only performed when \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/pollen-glia}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Pollen A et al. (2017). 
#' Molecular Identity of Human Outer Radial Glia during Cortical Development
#' \emph{Cell} 163, 55-67.
#'
#' @examples
#' sce <- PollenGliaData()
#' 
#' @export
PollenGliaData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.6.0"

    sce <- .create_sce(file.path("pollen-glia", version), has.rowdata=FALSE)

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
