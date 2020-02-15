#' Obtain the Zeisel brain data
#'
#' Obtain the mouse brain single-cell RNA-seq dataset from Zeisel et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' 
#' @details
#' Row data contains a single \code{"featureType"} field describing the type of each feature
#' (endogenous genes, mitochondrial genes, spike-in transcripts and repeats).
#' Spike-ins and repeats are stored as separate entries in the \code{\link{altExps}}.
#'
#' Column metadata is provided in the same form as supplied in \url{http://linnarssonlab.org/cortex/}.
#' This contains information such as the cell diameter and the published cell type annotations.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zeisel-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Zeisel A et al. (2015). 
#' Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. 
#' \emph{Science} 347(6226), 1138-42.
#'
#' @examples
#' sce <- ZeiselBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
#' @importFrom SummarizedExperiment rowData
ZeiselBrainData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("zeisel-brain", version))

    status <- rowData(sce)$featureType
    status[status=="mito"] <- "endogenous"
    sce <- splitAltExps(sce, status, "endogenous")

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Mm",
        ensembl=ensembl,
        location=location)
}
