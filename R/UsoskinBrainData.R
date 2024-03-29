#' Obtain the Usoskin brain data
#'
#' Obtain the mouse brain single-cell RNA-seq dataset from Usoskin et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided in the same form as supplied in External Table 2 of \url{http://linnarssonlab.org/drg/}.
#' This contains information such as the library of origin and the cell type.
#'
#' The count matrix contains information for repeats, marked with \code{r_} prefixes in the row names;
#' as well as mitochondrial transcripts, marked with \code{mt-} prefixes.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/usoskin-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of RPMs.
#'
#' @author Aaron Lun
#'
#' @references
#' Usoskin A et al. (2015).
#' Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. 
#' \emph{Nat. Neurosci.} 18(1), 145-53.
#'
#' @examples
#' sce <- UsoskinBrainData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData colData
UsoskinBrainData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("usoskin-brain-2015", "2023-12-19", realize.assays=TRUE)
    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("usoskin-brain", version), assays="rpm")
        rownames(sce) <- rowData(sce)[,1]
    }

    colnames(sce) <- colData(sce)[["Sample ID"]]
    
    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Mm",
        ensembl=ensembl,
        location=location)
}
