#' Obtain the Lawlor pancreas data
#'
#' Provides the human pancreas single-cell RNA-seq data from Lawlor et al. (2017).
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE86469. 
#' This contains information such as the cell type labels and patient status.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/lawlor-pancreas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Lawlor N et al. (2017). 
#' Single-cell transcriptomes identify human islet cell signatures and reveal cell-type-specific expression changes in type 2 diabetes.
#' \emph{Genome Res.} 27(2), 208-222.
#'
#' @examples
#' sce <- LawlorPancreasData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
LawlorPancreasData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("lawlor-pancreas", version), has.rowdata=FALSE)
    sce
}
