#' Obtain the Richard T cell data
#'
#' Obtain the mouse CD8+ T cell single-cell RNA-seq data from Richard et al. (2018).
#'
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in E-MTAB-6051.
#' This contains information such as the stimulus, time after stimulation, age of the mice and sequencing batch.
#'
#' Count data for ERCC spike-ins are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
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
RichardTCellData <- function(location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("richard-tcell-2018", "2023-12-19", realize.assays=TRUE)

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("richard-tcell", version), has.rowdata=FALSE)

        spike.type <- ifelse(grepl("ERCC", rownames(sce)), "ERCC", "endogenous")
        sce <- splitAltExps(sce, spike.type, ref="endogenous")

        spike.exp <- altExp(sce, "ERCC")
        spikedata <- ERCCSpikeInConcentrations(volume = 1000, dilution = 3e07)
        spikedata <- spikedata[rownames(spike.exp), ]
        rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
        altExp(sce, "ERCC") <- spike.exp
    }

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
