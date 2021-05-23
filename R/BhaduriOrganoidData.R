#' Obtain the Bhaduri cortical organoid data
#'
#' Obtain the human cortical organoid single-cell RNA-seq dataset from Bhaduri et al. (2020).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column data contains sample-level information.
#' In theory, there is also cell-level metadata for this dataset but it could not be unambiguously mapped to the column names.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bhaduri-organoid}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of normalized expression values.
#'
#' @author Aaron Lun
#'
#' @references
#' Bhaduri A et al. (2020). 
#' Cell stress in cortical organoids impairs molecular subtype specification.
#' \emph{Nature} 578(7793), 142-148.
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- BhaduriOrganoidData()
#' }
#' @export
#' @importFrom SingleCellExperiment splitAltExps altExp altExp<- reducedDims<- colPairs<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics cbind
BhaduriOrganoidData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.6.0"
    sce <- .create_sce(file.path("bhaduri-organoid", version), has.rowdata=FALSE)
    assayNames(sce) <- "normalized"

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
