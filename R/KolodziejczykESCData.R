#' Obtain the Kolodziejcyzk ESC data
#'
#' Obtain the mouse embryonic stem cell single-cell RNA-seq data from Kolodziejczyk et al. (2015).
#'
#' @param remove.htseq Logical scalar indicating whether HT-seq alignment statistics should be removed.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is generated from the column names,
#' and contains the culture conditions and the plate of origin for each cell.
#'
#' Count data for ERCC spike-ins are stored in the \code{"ERCC"} entry in the \code{\link{altExps}}.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/kolodziejczyk-esc}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun 
#'
#' @references
#' Messmer T et al. (2019). 
#' Transcriptional heterogeneity in naive and primed human pluripotent stem cells at single-cell resolution. 
#' \emph{Cell Rep} 26(4), 815-824.e4
#'
#' @examples
#' sce <- KolodziejczykESCData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
KolodziejczykESCData <- function(remove.htseq=TRUE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("kolodziejczyk-esc-2015", "2023-12-17", realize.assays=TRUE)

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("kolodziejczyk-esc", version), has.rowdata=FALSE, has.coldata=FALSE)

        sce$culture <- sub(".*mES_([^_]+)_.*", "\\1", colnames(sce))
        sce$plate <- sub(".*mES_[^_]+_([^_]+)_.*", "\\1", colnames(sce))

        if (remove.htseq) {
            sce <- sce[grep("^__", rownames(sce), invert=TRUE),]
        }

        spike.type <- ifelse(grepl("ERCC", rownames(sce)), "ERCC", "endogenous")
        sce <- splitAltExps(sce, spike.type, ref="endogenous")

        spike.exp <- altExp(sce, "ERCC")
        spikedata <- ERCCSpikeInConcentrations(volume = 1, dilution = 25e6)
        spikedata <- spikedata[rownames(spike.exp), ]
        rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
        altExp(sce, "ERCC") <- spike.exp
    }

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
