#' Obtain the Richard T cell data
#'
#' Obtain the mouse CD8+ T cell single-cell RNA-seq data from Richard et al. (2018).
#'
#' @details
#' Column metadata is provided in the same form as supplied in E-MTAB-6051.
#' This contains information such as the stimulus, time after stimulation, age of the mice and sequencing batch.
#'
#' Count data for ERCC spike-ins are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/richard-tcell}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Richard AC et al. (2018). 
#' T cell cytolytic capacity is independent of initial stimulation strength. 
#' \emph{Nat. Immunol.} 19(8), 849-858.
#'
#' @examples
#' sce <- RichardTCellData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
RichardTCellData <- function(location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("richard-tcell", version), has.rowdata=FALSE)

    spike.type <- ifelse(grepl("ERCC", rownames(sce)), "ERCC", "endogenous")
    sce <- splitAltExps(sce, spike.type, ref="endogenous")

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
