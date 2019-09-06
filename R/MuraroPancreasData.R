#' Obtain the Muraro pancreas data
#'
#' Download and cache the Muraro pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @param ensembl Logical scalar indicating whether the row names of the returned object should contain Ensembl identifiers.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#' 
#' @details
#' This function provides the human pancreas scRNA-seq data from Muraro et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts. 
#'
#' Row data contains fields for the symbol and chromosomal location of each gene.
#'
#' Column metadata is derived from the columns of the count matrix provided in GSE85241,
#' with additional cell type labels obtained from the authors (indirectly, via the Hemberg group).
#' Some cells have \code{NA} labels and were presumably removed prior to downstream analyses.
#'
#' ERCC spike-ins are represented as an alternative experiment.
#'
#' @return A \linkS4class{SingleCellExperiment} object.
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
        sce <- .convert_to_ensembl(sce, rowData(sce)$symbol, species="Hs")
    }
    sce
}
