#' Obtain the Xin pancreas data
#'
#' Obtain the human pancreas single-cell RNA-seq dataset from Xin et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Row data contains fields for the Entrez ID and symbol for each gene.
#' Column metadata was obtained from the authors (indirectly, via the Hemberg group) 
#' and contains information such as the cell type labels and donor status.
#'
#' If \code{ensembl=TRUE}, the Entrez IDs are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/xin-pancreas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of RPKMs.
#'
#' @author Aaron Lun,
#' using additional metadata obtained by Vladimir Kiselev.
#'
#' @references
#' Xin A et al. (2016).
#' RNA sequencing of single human islet cells reveals type 2 diabetes genes.
#' \emph{Cell Metab.} 24(4), 608-615.
#'
#' @examples
#' sce <- XinPancreasData()
#' 
#' @export
XinPancreasData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("xin-pancreas", version), assays="rpkm")

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs", 
        keytype="ENTREZID",
        ensembl=ensembl,
        location=location)
}
