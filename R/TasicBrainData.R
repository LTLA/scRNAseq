#' Obtain the Tasic brain data
#'
#' Obtain the mouse brain single-cell RNA-seq data from Tasic et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata is provided in the same form as supplied in GSE71585.
#' This contains information such as the reporter gene expressed in each cell, the mouse line, dissection type and so on.
#'
#' Count data for ERCC spike-ins are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#' Note that some of the spike-in rows have \code{NA} observations for some (but not all) cells.
#'
#' The last 9 columns (containing \code{_CTX_} in their names) correspond to no-cell control libraries.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/tasic-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Tasic A et al. (2016). 
#' Adult mouse cortical cell taxonomy revealed by single cell transcriptomics.
#' \emph{Nat. Neurosci.} 19(2), 335-46.
#'
#' @examples
#' sce <- TasicBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
TasicBrainData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("tasic-brain", version), has.rowdata=FALSE)

    status <- ifelse(grepl("^ERCC-[0-9]+$", rownames(sce)), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")
    spike.exp <- altExp(sce, "ERCC")
    spikedata <- ERCCSpikeInConcentrations(volume = 100, dilution = 1000000)
    spikedata <- spikedata[rownames(spike.exp), ]
    rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
    altExp(sce, "ERCC") <- spike.exp

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Mm",
        ensembl=ensembl,
        location=location)
}
