#' Obtain the Chen brain data
#'
#' Obtain the mouse brain single-cell RNA-seq data from Chen et al. (2017).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE87544.
#' This contains the putative cell type assigned by the original authors.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/chen-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Chen R et al. (2017).
#' Single-Cell RNA-Seq reveals hypothalamic cell diversity. 
#' \emph{Cell Rep.} 18, 3227-3241.
#'
#' @examples
#' sce <- ChenBrainData()
#' 
#' @export
ChenBrainData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("chen-brain-2017", "2023-12-14", realize.assays=TRUE)
    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("chen-brain", version), has.rowdata=FALSE)
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
