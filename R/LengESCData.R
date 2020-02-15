#' Obtain the Leng ESC data
#'
#' Obtain the human embryonic stem cell single-cell RNA-seq data from Leng et al. (2015).
#'
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Column metadata contains the cell line, experiment number and experimentally determined cell cycle phase for each cell,
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/leng-esc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of normalized expected read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Leng F et al. (2015).
#' Oscope identifies oscillatory genes in unsynchronized single-cell RNA-seq experiments. 
#' \emph{Nat. Methods} 12(10), 947-950.
#' 
#' @examples
#' sce <- LengESCData()
#' 
#' @export
LengESCData <- function(location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("leng-esc", version), assays="normcounts", has.rowdata=FALSE)
    .define_location_from_ensembl(sce, species="Hs", location=location)
}
