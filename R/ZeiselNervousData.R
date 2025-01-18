#' Obtain the Zeisel nervous system data
#'
#' Obtain the mouse nervous system single-cell RNA-seq dataset from Zeisel et al. (2018).
#'
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' Row data contains the gene symbol as well as some relevant per-gene statistics,
#' e.g., the squared coefficient of variance, mean, and whether it was selected for downstream analyses.
#'
#' Column data contains a wide variety of fields including patient-level information, sample-level sequencing statistics and many flavors of cell type classification.
#' Note that many numeric columns may have \code{NA} values if they could not be successfully parsed form the source file.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link[SummarizedExperiment]{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zeisel-nervous}.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Zeisel A et al. (2018). 
#' Molecular architecture of the mouse nervous system.
#' \emph{Cell} 174(4), 999-1014.
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- ZeiselNervousData()
#' }
#' @export
#' @importFrom SingleCellExperiment splitAltExps altExp altExp<- reducedDims<- colPairs<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics cbind
ZeiselNervousData <- function(location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("zeisel-nervous-2018", "2023-12-22", realize.assays=TRUE)

    } else {
        version <- "2.6.0"

        hub <- ExperimentHub()
        sce <- .create_sce(file.path("zeisel-nervous", version), hub=hub)

        reducedDims(sce) <- hub[[which(hub$rdatapath==file.path("scRNAseq/zeisel-nervous", version, "reddims.rds"))]]
        colPairs(sce) <- hub[[which(hub$rdatapath==file.path("scRNAseq/zeisel-nervous", version, "colpairs.rds"))]]
    }

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
