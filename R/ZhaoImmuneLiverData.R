#' Obtain the Zhao immune liver data
#'
#' Obtain the human liver immune single-cell RNA-seq data from Zhao et al. (2020).
#'
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param filter Logical scalar indicating if the filtered subset should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata contains various cell labels as provided by the authors.
#' Some of these labels required assembly on our part:
#' \itemize{
#' \item The \code{broad} label was assigned to each barcode based on whether that barcode was present in each \code{*_identities.tsv.gz} in GSE125188's supplementary files.
#' \item For each cell barcode that was present in one of these files, the \code{fine} label was generated from the \code{Group} annotations inside that file.
#' }
#' We guessed the \code{sample} for each cell by assuming that the GEM group numbers match the order of samples in GSE125188.
#' We also assumed that \dQuote{donor 4} is a typo, given that the paper only mentions 3 donors.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link[SummarizedExperiment]{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#' 
#' If \code{filter=TRUE}, only cells that have been used in the original analysis are returned.
#' Otherwise, the cells used are specified in the \code{retained} column of the \code{\link[SummarizedExperiment]{colData}}.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zhao-immune-liver}.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Zhao J et al. (2020). 
#' Single-cell RNA sequencing reveals the heterogeneity of liver-resident immune cells in human.
#' \emph{Cell Discov} 6, 22.
#'
#' @examples
#' sce.zhao <- ZhaoImmuneLiverData()
#' 
#' @export
ZhaoImmuneLiverData <- function(location=TRUE, filter=FALSE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("zhao-immune-2020", "2023-12-22", realize.assays=TRUE)

    } else {
        version <- "2.6.0"
        sce <- .create_sce(file.path("zhao-immune-liver", version), has.rowdata=TRUE)

        if (filter) {
            sce <- sce[,sce$retained]
            sce$retained <- NULL
        }

        gem.group <- as.integer(sub(".*-", "", colnames(sce)))
        sce$sample <- c("donor 1 blood", "donor 1 spleen", "donor 1 liver",
            "donor 2 blood", "donor 2 spleen", "donor 2 liver",
            "donor 3 blood", "donor 3 spleen", "donor 3 liver")[gem.group]
    }

    sce
}
