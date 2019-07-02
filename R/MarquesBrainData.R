#' Obtain the Marques brain data
#'
#' Download and cache the Marques brain single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the brain scRNA-seq data from Marques et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' Column metadata is provided in the same form as supplied in GSE75330. 
#' This contains information such as the cell type and age/sex of the mouse of origin for each cell.
#' 
#' Note that some genes may be present in multiple rows corresponding to different genomic locations.
#' These additional rows are identified by a \code{_loc[2-9]} prefix in their row names.
#' Users may wish to consider either removing them or merging them, e.g., with \code{scater::sumCountsAcrossFeatures}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Marques A et al. (2016).
#' Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system. 
#' \emph{Science} 352(6291), 1326-1329. 
#'
#' @examples
#' sce <- MarquesBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike
MarquesBrainData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("tasic-brain", version), assays="counts", has.rowdata=FALSE)
    isSpike(sce, "ERCC") <- grep("^ERCC-[0-9]+", rownames(sce))
    sce
}
