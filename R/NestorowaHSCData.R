#' Obtain the Nestorowa HSC data
#'
#' Obtain the mouse haematopoietic stem cell single-cell RNA-seq data from Nestorowa et al. (2015).
#'
#' @param remove.htseq Logical scalar indicating whether HT-seq alignment statistics should be removed.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Rows corresponding to HT-seq's alignment statistics are removed by default.
#' These can be retained by setting \code{remove.htseq=FALSE}.
#'
#' Column metadata includes the cell type mapping, as described on the website (see References),
#' and the FACS expression levels of selected markers.
#' Note that these are stored as nested matrices within the \code{\link{colData}}.
#'
#' Diffusion map components are provided as the \code{"diffusion"} entry in the \code{\link{reducedDims}}.
#'
#' Counts for ERCC spike-ins are stored in the \code{"ERCC"} entry in the \code{\link{altExps}}.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/nestorowa-hsc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
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
#' @importFrom SingleCellExperiment splitAltExps reducedDim<-
NestorowaHSCData <- function(remove.htseq=TRUE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("nestorowa-hsc", version), has.rowdata=FALSE)

    if (remove.htseq) {
        sce <- sce[grep("^__", rownames(sce), invert=TRUE),]
    }

    reducedDim(sce, "diffusion") <- sce$diffusion
    sce$diffusion <- NULL

    status <- ifelse(grepl("^ERCC-[0-9]+", rownames(sce)), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
