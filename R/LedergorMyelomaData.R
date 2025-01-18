#' Obtain the Ledergor Myeloma data
#'
#' Obtain the human multiple myeloma single-cell RNA-seq data from Ledergor et al. (2018).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata was created from the sample metadata file in GSE117156. It contains an 'Experiment_ID' column,
#' from which the tissue and subject of origin were extracted, as well as the condition and treatment status of the subject.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link[SummarizedExperiment]{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ledergor-myeloma}.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Milan Malfait
#'
#' @references
#' Ledergor G et al. (2018)
#' Single cell dissection of plasma cell heterogeneity in symptomatic and
#' asymptomatic myeloma. \emph{Nat Med.} 24(12), 1867–1876.
#'
#' @examples
#' sce <- LedergorMyelomaData()
#'
#' @export
#' @importFrom SummarizedExperiment rowData
LedergorMyelomaData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("ledergor-myeloma-2018", "2023-12-20", realize.assays=TRUE)
    } else {
        version <- "2.8.0"
        sce <- .create_sce(file.path("ledergor-myeloma", version), has.rowdata=FALSE)
    }

    .convert_to_ensembl(sce,
        species = "Hs",
        symbols = rownames(sce),
        ensembl = ensembl,
        location = location)
}
