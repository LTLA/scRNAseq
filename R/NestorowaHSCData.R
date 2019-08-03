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
#' Row corresponding to HT-seq's alignment statistics are present in the count matrix.
#' By default, these are removed prior to downstream analyses.
#'
#' Column metadata includes the cell type mapping, as described on the website (see References),
#' the diffusion map components and the FACS expression levels of selected markers.
#'
#' ERCC spike-ins are stored as alternative experiments.
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
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment splitSCEByAlt
NestorowaHSCData <- function(remove.htseq=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("nestorowa-hsc", version), has.rowdata=FALSE)

    if (remove.htseq) {
        sce <- sce[grep("^__", rownames(sce), invert=TRUE),]
    }

    status <- ifelse(grepl("^ERCC-[0-9]+", rownames(sce)), "ERCC", "endogenous")
    splitSCEByAlt(sce, status, ref="endogenous")
}
