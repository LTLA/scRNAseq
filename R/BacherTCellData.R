#' Obtain the Bacher T cell data
#'
#' Obtain the human COVID T cell single-cell RNA-seq dataset from Bacher et al. (2020).
#'
#' @param filtered Logical scalar indicating whether to filter out cells that were not used by the authors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' Column metadata is scraped from GEO, using both the author-supplied TSV of per-cell annotations and the sample-level metadata.
#' This contains information such as the diagnosis, severity, WHO class, clustering and clonotype.
#'
#' If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
#' Otherwise, an additional \code{filtered} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors. 
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bacher-tcell}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Bacher P et al. (2020). 
#' Low avidity T cell responses to SARS-CoV-2 in unexposed individuals and severe COVID-19
#' \emph{Immunity} 53, 1258-1271
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- BacherTCellData()
#' }
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps altExp altExp<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics cbind
BacherTCellData <- function(filtered=TRUE, ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("bacher-tcell-2020", "2023-12-21", realize.assays=TRUE)
    } else {
        version <- "2.6.0"
        sce <- .create_sce(file.path("bacher-tcell", version), has.rowdata=FALSE)
    }

    if (filtered) {
        sce <- sce[,sce$retained]
        sce$retained <- NULL
    }

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
