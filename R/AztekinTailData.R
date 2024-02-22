#' Obtain the Aztekin tail data
#'
#' Obtain the Xenopus tail single-cell RNA-seq data from Aztekin et al. (2019).
#'
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in E-MTAB-7761.
#' This contains information such as the treatment condition, batch, putative cell type, putative cell cycle phase.
#'
#' The UMAP results are available as the \code{"UMAP"} entry in the \code{\link{reducedDims}}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/aztekin-tail}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Aztekin C et al. (2019).
#' Identification of a regeneration-organizing cell in the Xenopus tail.
#' \emph{Science} 364(6441), 653-658
#'
#' @examples
#' sce <- AztekinTailData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps reducedDim<- 
#' @importFrom SummarizedExperiment colData colData<-
AztekinTailData <- function(legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("aztekin-tail-2019", "2023-12-14", realize.assays=TRUE)

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("aztekin-tail", version), has.rowdata=FALSE)

        # Move UMAP coordinates to the reducedDims.
        reducedDim(sce, "UMAP") <- as.matrix(colData(sce)[,c("X", "Y")])
        colData(sce) <- colData(sce)[,setdiff(colnames(colData(sce)), c("X", "Y"))]
    }

    sce
}
