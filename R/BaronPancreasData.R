#' Obtain the Baron pancreas data
#'
#' Obtain the human/mouse pancreas single-cell RNA-seq data from Baron et al. (2017).
#'
#' @param which String specifying the species to get data for.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE84133.
#' This contains information such as the cell type labels and donor ID (for humans) or strain (for mouse).
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/baron-pancreas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Baron M et al. (2017). 
#' Single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. 
#' \emph{Cell Syst.} 3(4), 346-360.
#'
#' @examples
#' sce.human <- BaronPancreasData()
#'
#' sce.mouse <- BaronPancreasData("mouse")
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
BaronPancreasData <- function(which=c("human", "mouse"), ensembl=FALSE, location=TRUE, legacy=FALSE) {
    which <- match.arg(which)

    if (!legacy) {
        sce <- fetchDataset("baron-pancreas-2016", which, realize.assays=TRUE)
    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("baron-pancreas", version), has.rowdata=FALSE, suffix=which)
    }

    .convert_to_ensembl(sce, 
        ensembl=ensembl,
        species=if (which=="human") "Hs" else "Mm",
        symbols=rownames(sce),
        location=location)
}
