#' Obtain the Buettner ESC data
#'
#' Obtain the mouse embryonic stem cell single-cell RNA-seq data from Buettner et al. (2015).
#'
#' @param remove.htseq Logical scalar indicating whether HT-seq alignment statistics should be removed.
#'
#' @details
#' Rows corresponding to HT-seq's alignment statistics are removed by default.
#' These can be retained by setting \code{remove.htseq=FALSE}.
#'
#' Column metadata contains the experimentally determined cell cycle phase for each cell.
#'
#' Counts for ERCC spike-ins are stored in the \code{"ERCC"} entry in the \code{\link{altExps}}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/buettner-esc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Buettner F et al. (2015).
#' Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells.
#' \emph{Nat. Biotechnol.} 33(2), 155-160.
#' 
#' @examples
#' sce <- BuettnerESCData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps 
BuettnerESCData <- function(remove.htseq=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("buettner-esc", version), has.rowdata=TRUE)

    if (remove.htseq) {
        sce <- sce[grepl("^ENSMUS", rownames(sce)) | grepl("^ERCC", rownames(sce)),]
    }

    status <- ifelse(grepl("^ERCC-[0-9]+", rownames(sce)), "ERCC", "endogenous")
    splitAltExps(sce, status, ref="endogenous")
}
