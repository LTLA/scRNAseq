#' Obtain the Marques brain data
#'
#' Obtain the mouse brain single-cell RNA-seq data from Marques et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE75330. 
#' This contains information such as the cell type and age/sex of the mouse of origin for each cell.
#' 
#' Note that some genes may be present in multiple rows corresponding to different genomic locations.
#' These additional rows are identified by a \code{_loc[2-9]} suffix in their row names.
#' Users may wish to consider either removing them or merging them, e.g., with \code{scater::sumCountsAcrossFeatures}.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#' All searching is performed after removing the \code{_loc[2-9]} suffix. 
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/marques-brain}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
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
MarquesBrainData <- function(ensembl=FALSE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("marques-brain", version), assays="counts", has.rowdata=FALSE)

    if (ensembl) {
        sce <- .convert_to_ensembl(sce, species="Mm", symbols=sub("_loc[0-9]+$", "", rownames(sce)))
    }
    sce
}
