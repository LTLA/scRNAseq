#' Obtain the Grun pancreas data
#'
#' Download and cache the Grun pancreas single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @details
#' This function provides the human pancreas scRNA-seq data from Grun et al. (2016)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts. 
#'
#' Row data contains fields for the symbol and chromosomal location of each gene.
#' Spike-ins are specially labelled with the \code{\link{isSpike}} function.
#'
#' Column metadata is derived from the columns of the count matrix provided in GSE85241,
#' with additional cell type labels obtained from the authors (indirectly, via the Hemberg group) 
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun,
#' using additional metadata obtained by Vladimir Kiselev.
#'
#' @references
#' Grun D et al. (2016). 
#' De novo prediction of stem cell identity using single-cell transcriptome data. 
#' \emph{Cell Stem Cell} 19(2), 266-277. 
#'
#' @examples
#' sce <- GrunPancreasData()
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment isSpike<-
GrunPancreasData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("grun-pancreas", version), has.rowdata=FALSE) 
    symbol <- sub("__.*", "", rownames(sce))
    loc <- sub(".*__", "", rownames(sce))
    rowData(sce) <- DataFrame(symbol=symbol, chr=loc)
    isSpike(sce, "ERCC") <- grep("ERCC-[0-9]+", symbol)
    sce
}
