#' Obtain the Ernst spermatogenesis data
#'
#' Obtain the mouse spermatogenesis single-cell RNA-seq dataset from Ernst et al. (2019).
#'
#' @param method String indicating which cell caller to obtain results for.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' This study contains two analyses done with datasets from different cell calling algorithms.
#' One uses Cellranger version 2 while the other uses \code{emptyDrops} from \pkg{DropletUtils}.
#'
#' Column metadata includes sample information, per-cell QC metrics and cell type labels.
#' In particular, the sample label specifies the developmental stage of the mouse.
#'
#' Note that \code{method="Cellranger"} contains additional data for Tc1 mice.
#' These mice have an additional human chromosome 21 inserted alongside the usual mouse chromosomes.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ernst-spermatogenesis}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Ernst C et al. (2019). 
#' Staged developmental mapping and X chromosome transcriptional dynamics during mouse spermatogenesis.
#' \emph{Nat Commun} 10, 1251
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- ErnstSpermatogenesisData()
#' }
#' 
#' @export
ErnstSpermatogenesisData <- function(method=c("emptyDrops", "Cellranger"), location=TRUE) {
    version <- "2.6.0"

    method <- match.arg(method)
    sce <- .create_sce(file.path("ernst-spermatogenesis", version), suffix=tolower(method))

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
