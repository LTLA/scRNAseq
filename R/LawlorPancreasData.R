#' Obtain the Lawlor pancreas data
#'
#' Download and cache the Lawlor pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the pancreas scRNA-seq data from Lawlor et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Column metadata is provided in the same form as supplied in GSE86469. 
#' This contains information such as the cell type labels and patient status.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
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
