#' Obtain the Nestorowa HSC data
#'
#' Download and cache the Nestorowa HSC single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @param remove.htseq Logical scalar indicating whether HT-seq alignment statistics should be removed.
#'
#' @details
#' This function provides the haematopoietic stem cell scRNA-seq data from Nestorowa et al. (2015)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Rows containing spike-in transcripts are specially labelled with the \code{\link{isSpike}} function.
#' Note that rows corresponding to HT-seq's alignment statistics are also included in the design matrix;
#' these should probably be removed prior to downstream analyses.
#'
#' Column metadata includes the cell type mapping, as described on the website (see References),
#' the diffusion map components and the FACS expression levels of selected markers.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Nestorowa S et al. (2016). 
#' A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation
#' \emph{Blood} 128, e20-e31.
#'
#' Gene and protein expression in adult haematopoiesis: Data.
#' \url{http://blood.stemcells.cam.ac.uk/single_cell_atlas.html#data}.
#' 
#' @examples
#' sce <- NestorowaHSCData()
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
#' @importFrom SummarizedExperiment rowData
NestorowaHSCData <- function(remove.htseq=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("nestorowa-hsc", version), has.rowdata=FALSE)
    if (remove.htseq) {
        sce <- sce[grep("^__", rownames(sce), invert=TRUE),]
    }
    isSpike(sce, "ERCC") <- grep("^ERCC-[0-9]+", rowData(sce)$refseq)
    sce
}
