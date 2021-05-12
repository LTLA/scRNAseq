#' Obtain the Jessa brain data
#'
#' Obtain the mouse brain single-cell RNA-seq dataset from Jessa et al. (2019).
#'
#' @param filtered Logical scalar indicating whether to filter out cells that were not used by the authors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
#' Otherwise, an additional \code{filtered} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors. 
#'
#' The column data contains sample of origin, some QC metrics and various cluster assignments for each cell.
#' Cluster assignments starting with \code{Sample_*} are derived from per-sample analyses and cannot be compared sensibly across samples.
#' Other clusterings (\code{Forebrain_*} and \code{Pons_*}) are derived from joint analyses across all samples involving the named tissue.
#'
#' The \code{reducedDims} of the output contains various dimensionality reduction results.
#' Coordinates for entries prefixed with \code{Sample_*} were generated from per-sample analyses and cannot be compared across samples.
#' Coordinates for entries prefixed with \code{Forebrain_*} and \code{Pons_*} were generated from joint analyses from the corresponding tissue. 
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/jessa-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Jessa S et al. (2019).
#' Stalled developmental programs at the root of pediatric brain tumors
#' \emph{Nat Genet} 51, 1702-1713
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- JessaBrainData()
#' } 
#' @export
#' @importFrom SingleCellExperiment reducedDims<-
JessaBrainData <- function(filtered=TRUE, location=TRUE) {
    version <- "2.6.0"

    hub <- .ExperimentHub()
    sce <- .create_sce(file.path("jessa-brain", version), hub=hub)

    reddims <- hub[hub$rdatapath==file.path("scRNAseq", "jessa-brain", version, "reddims.rds")][[1]]
    reducedDims(sce) <- reddims

    if (filtered) {
        sce <- sce[,sce$retained]
        sce$filtered <- NULL
    }

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
