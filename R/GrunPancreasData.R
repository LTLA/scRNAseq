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
#' Column metadata is derived from the column names of the count matrix with the sample annotations in GSE81076.
#' This includes the donor identity for each cell and the type of sample.
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
#' @importFrom SummarizedExperiment rowData<- colData<-
#' @importFrom SingleCellExperiment isSpike<-
GrunPancreasData <- function() {
    version <- "2.0.0"
    sce <- .create_sce(file.path("grun-pancreas", version), has.rowdata=FALSE, has.coldata=FALSE)

    # Loading the column metadata
    lib.names <- sub("_.*", "", colnames(counts))
    donor.names <- sub("(D10|D17|D2|D3|D7).*", "\\1", lib.names)
    treatment <- c(
        D10631="CD63+ sorted cells",
        D101="live sorted cells, library 1",
        D102="live sorted cells, library 2",
        D1713="CD13+ sorted cells",
        D172444="CD24+ CD44+ live sorted cells",
        D17All1="live sorted cells, library 1",
        D17All2="live sorted cells, library 2",
        D17TGFB="TGFBR3+ sorted cells",
        D2ex="exocrine fraction, live sorted cells",
        D3en1="live sorted cells, library 1",
        D3en2="live sorted cells, library 2",
        D3en3="live sorted cells, library 3",
        D3en4="live sorted cells, library 4",
        D3ex="exocrine fraction, live sorted cells",
        D71="live sorted cells, library 1",
        D72="live sorted cells, library 2",
        D73="live sorted cells, library 3",
        D74="live sorted cells, library 4"
    )[lib.names]
    colData(sce) <- DataFrame(donor=donor.names, sample=treatment)

    # Splitting up gene information.
    symbol <- sub("__.*", "", rownames(sce))
    loc <- sub(".*__", "", rownames(sce))
    rowData(sce) <- DataFrame(symbol=symbol, chr=loc)
    isSpike(sce, "ERCC") <- grep("ERCC-[0-9]+", symbol)
    sce
}
