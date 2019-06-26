#' Obtain the La Manno brain data
#'
#' Download and cache the La Manno brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @param which A string specifying which data set should be obtained.
#' 
#' @details
#' This function provides the brain scRNA-seq data from La Manno et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with one of several possible matrices of UMI counts.
#'
#' Column metadata is provided in the same form as supplied in the supplementary tables in GSE71585.
#' This contains information such as the time point and cell type.
#'
#' The various settings of \code{which} will obtain different data sets.
#' \itemize{
#' \item \code{"human-es"}, human embryonic stem cells.
#' \item \code{"human-embryo"}, human embryo midbrain.
#' \item \code{"human-ips"}, human induced pluripotent stem cells.
#' \item \code{"mouse-adult"}, mouse adult dopaminergic neurons.
#' \item \code{"mouse-embryo"}, mouse embryo midbrain.
#' }
#' Unfortunately, each of these datasets uses a different set of features.
#' If multiple datasets are to be used simultaneously, users will have to decide how to merge them,
#' e.g., by taking the intersection of common features across all datasets.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' La Manno A et al. (2016).
#' Molecular diversity of midbrain development in mouse, human, and stem cells. 
#' \emph{Cell} 167(2), 566-580.
#'
#' @examples
#' sce <- LaMannoBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike<-
LaMannoBrainData <- function(which=c("human-es", "human-embryo", "human-ips", "mouse-adult", "mouse-embryo")) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("lamanno-brain", version), has.rowdata=FALSE, suffix=which)
    sce
}
