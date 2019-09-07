#' Obtain the Muraro pancreas data
#'
#' Obtain the human pancreas single-cell RNA-seq data from Muraro et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' 
#' @details
#' Row data contains fields for the symbol and chromosomal location of each gene.
#'
#' Column metadata is derived from the columns of the count matrix provided in GSE85241,
#' with additional cell type labels obtained from the authors (indirectly, via the Hemberg group).
#' Some cells have \code{NA} labels and were presumably removed prior to downstream analyses.
#'
#' ERCC spike-ins are represented as an alternative experiment.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/muraro-pancreas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun,
#' using additional metadata obtained by Vladimir Kiselev.
#'
#' @references
#' Muraro MJ et al. (2016). 
#' A single-cell transcriptome atlas of the human pancreas.
#' \emph{Cell Syst.} 3(4), 385-394.
#'
#' @examples
#' sce <- MuraroPancreasData()
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment splitAltExps
MuraroPancreasData <- function(ensembl=FALSE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("muraro-pancreas", version), has.rowdata=FALSE) 

    symbol <- sub("__.*", "", rownames(sce))
    loc <- sub(".*__", "", rownames(sce))
    rowData(sce) <- DataFrame(symbol=symbol, chr=loc)

    status <- ifelse(grepl("^ERCC-[0-9]+", symbol), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")

    if (ensembl) {
        sce <- .convert_to_ensembl(sce, symbols=rowData(sce)$symbol, species="Hs")
    }
    sce
}
