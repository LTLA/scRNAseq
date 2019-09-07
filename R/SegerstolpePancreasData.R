#' Obtain the Segerstolpe pancreas data
#'
#' Download the human pancreas single-cell RNA-seq (scRNA-seq) dataset from Segerstolpe et al. (2016)
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#'
#' @details
#' Row data contains fields for the gene symbol and RefSeq transcript IDs corresponding to each gene.
#' The rows of the output object are named with the symbol, but note that these are not unique.
#'
#' Column metadata were extracted from the \code{Characteristics} fields of the SDRF file for ArrayExpress E-MTAB-5061.
#' This contains information such as the cell type labels and patient status.
#'
#' ERCC spike-ins are represented as an alternative experiment.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/segerstolpe-pancreas}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Segerstolpe A et al. (2016). 
#' Single-cell transcriptome profiling of human pancreatic islets in health and type 2 diabetes. 
#' \emph{Cell Metab.} 24(4), 593-607.
#'
#' @examples
#' sce <- SegerstolpePancreasData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
#' @importFrom SummarizedExperiment rowData
SegerstolpePancreasData <- function(ensembl=FALSE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("segerstolpe-pancreas", version))
    rownames(sce) <- rowData(sce)$symbol

    status <- ifelse(grepl("^ERCC-[0-9]+", rowData(sce)$refseq), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")

    if (ensembl) {
        sce <- .convert_to_ensembl(sce, symbols=rownames(sce), species="Hs")
    }
    sce
}
