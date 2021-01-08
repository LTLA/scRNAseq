#' Obtain the Darmanis brain data
#'
#' Obtain the human brain single-cell RNA-seq dataset from Darmanis et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata is scraped from GEO and includes patient information, tissue of origin and likely cell type. 
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' This is only performed when \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/bacher-tcell}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Darmanis S et al. (2020). 
#' A survey of human brain transcriptome diversity at the single cell level. 
#' \emph{Proc Natl Acad Sci USA} 112, 7285-90.
#'
#' @examples
#' sce <- DarmanisBrainData()
#' 
#' @export
DarmanisBrainData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.6.0"

    sce <- .create_sce(file.path("darmanis-brain", version), has.rowdata=FALSE)

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
