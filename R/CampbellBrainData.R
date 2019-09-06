#' Obtain the Campbell brain data
#'
#' Obtain the mouse brain single-cell RNA-seq data from Campbell et al. (2017).
#'
#' @param ensembl Logical scalar indicating whether the row names of the returned object should contain Ensembl identifiers.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE93374.
#' This contains information such as the diet of the mice, sex and proposed cell type for each cell.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/campbell-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Campbell R et al. (2017).
#' A molecular census of arcuate hypothalamus and median eminence cell types. 
#' \emph{Nat. Neurosci.} 20, 484-496. 
#'
#' @examples
#' sce <- CampbellBrainData()
#' 
#' @export
CampbellBrainData <- function(ensembl=FALSE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("campbell-brain", version), has.rowdata=FALSE)

    if (ensembl) {
        sce <- .convert_to_ensembl(sce, species="Mm", symbols=rownames(sce))
    }
    sce 
}
