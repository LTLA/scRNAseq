#' Obtain the Grun pancreas data
#'
#' Obtain the human pancreas single-cell RNA-seq data from Grun et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Row metadata contains fields for the symbol and chromosomal location of each gene,
#' as derived from the row names.
#'
#' Column metadata is derived from the column names of the count matrix with the sample annotations in GSE81076.
#' This includes the donor identity for each cell and the type of sample.
#'
#' The \code{"ERCC"} entry in the \code{\link{altExps}} contains count data for the ERCC spike-in transcripts.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/grun-pancreas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts. 
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
#' @importFrom SingleCellExperiment splitAltExps
GrunPancreasData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("grun-pancreas", version), has.rowdata=FALSE, has.coldata=FALSE)

    # Loading the column metadata
    lib.names <- sub("_.*", "", colnames(sce))
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
    colData(sce) <- DataFrame(donor=donor.names, sample=treatment, row.names=colnames(sce))

    # Splitting up gene information.
    symbol <- sub("__.*", "", rownames(sce))
    loc <- sub(".*__", "", rownames(sce))
    rowData(sce) <- DataFrame(symbol=symbol, chr=loc)

    # Splitting spike-ins into an alternative experiment.
    status <- ifelse(grepl("ERCC-[0-9]+", symbol), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")

    .convert_to_ensembl(sce, 
        symbols=rowData(sce)$symbol, 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
