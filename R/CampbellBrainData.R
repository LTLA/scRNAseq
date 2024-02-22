#' Obtain the Campbell brain data
#'
#' Obtain the mouse brain single-cell RNA-seq data from Campbell et al. (2017).
#'
#' @param ensembl Logical scalar indicating whether the row names of the returned object should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE93374.
#' This contains information such as the diet of the mice, sex and proposed cell type for each cell.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/campbell-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
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
CampbellBrainData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("campbell-brain-2017", "2023-12-14", realize.assays=TRUE)
    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("campbell-brain", version), has.rowdata=FALSE)
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
