#' Obtain the Bach mammary data
#'
#' Obtain the mouse mammary gland single-cell RNA-seq data from Bach et al. (2017).
#'
#' @param samples A character vector with at least one element, specifying which samples(s) to retrieve.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata is extracted from the sample annotation in GSE106273,
#' and refers to the developmental stage of the mammary gland.
#'
#' If multiple samples are specified in \code{samples}, the count matrices will be \code{cbind}ed together.
#' Cells originating from different samples are identifiable by the \code{"Sample"} field in the column metadata.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bach-mammary}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Bach K et al. (2017).
#' Differentiation dynamics of mammary epithelial cells revealed by single-cell RNA sequencing. 
#' \emph{Nat Commun.} 8(1), 2128
#'
#' @examples
#' sce <- BachMammaryData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom BiocGenerics cbind
BachMammaryData <- function(samples=c("NP_1", "NP_2", "G_1", "G_2", "L_1", "L_2", "PI_1", "PI_2"),
    location=TRUE)
{
    version <- "2.0.0"
    host <- file.path("bach-mammary", version)
    samples <- match.arg(samples, several.ok=TRUE)

    collected <- vector("list", length(samples))
    for (i in seq_along(samples)) {
       collected[[i]] <- .create_sce(host, has.rowdata=FALSE, suffix=samples[i])
    }
    sce <- do.call(cbind, collected)

    ehub <- ExperimentHub()
    rowData(sce) <- ehub[ehub$rdatapath==file.path("scRNAseq", host, "rowdata.rds")][[1]]

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
