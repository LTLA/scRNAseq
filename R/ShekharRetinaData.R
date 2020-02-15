#' Obtain the Shekhar retina data
#'
#' Obtain the mouse retina single-cell RNA-seq dataset from Shekhar et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata contains the cluster identities as reported in the paper.
#' Note that some cells will have \code{NA} identities as they are present in the count matrix but not in the metadata file.
#' These are presumably low-quality cells that were discarded prior to clustering.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/shekhar-retina}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Shekhar K et al. (2016). 
#' Comprehensive classification of retinal bipolar neurons by single-cell transcriptomics.
#' \emph{Cell} 166(5), 1308-1323.
#'
#' @examples
#' sce <- ShekharRetinaData()
#' 
#' @export
ShekharRetinaData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("shekhar-retina", version), has.rowdata=FALSE)

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
