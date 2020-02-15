#' Obtain the Paul HSC data
#'
#' Obtain the mouse haematopoietic stem cell single-cell RNA-seq data from Paul et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param discard.multiple Logical scalar indicating whether ambiguous rows should be discarded.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Column metadata includes the plate and the mouse of origin, fluoresence intensities from indexed sorting 
#' and the number of cells in each well.
#'
#' Some of the original rownames are concatenated symbols from multiple genes.
#' We consider these rows to represent ambiguously assigned counts and discard them if \code{discard.multiple=TRUE}.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/nestorowa-hsc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Paul F et al. (2015). 
#' Transcriptional heterogeneity and lineage commitment in myeloid progenitors. 
#' \emph{Cell} 163, 1663-77.
#'
#' @examples
#' sce <- PaulHSCData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
PaulHSCData <- function(ensembl=FALSE, discard.multiple=TRUE, location=TRUE) {
    version <- "2.2.0"
    sce <- .create_sce(file.path("paul-hsc", version), has.rowdata=FALSE)

    if (discard.multiple) {
        sce <- sce[!grepl(";", rownames(sce)),]
    }
    
    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
