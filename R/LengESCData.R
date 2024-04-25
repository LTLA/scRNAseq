#' Obtain the Leng ESC data
#'
#' Obtain the human embryonic stem cell single-cell RNA-seq data from Leng et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether gene symbols should be converted to Ensembl annotation.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata contains the cell line, experiment number and experimentally determined cell cycle phase for each cell.
#'
#' If \code{ensembl=TRUE}, the gene symbols in the published annotation are converted to Ensembl.
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
LengESCData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("leng-esc-2015", "2024-04-25", realize.assays=TRUE)
    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("leng-esc", version), assays="normcounts", has.rowdata=FALSE)
    }

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
