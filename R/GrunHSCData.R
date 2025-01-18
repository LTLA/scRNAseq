#' Obtain the Grun HSC data
#'
#' Obtain the mouse haematopoietic stem cell single-cell RNA-seq data from Grun et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Row metadata contains the symbol and chromosomal location for each gene.
#' Column metadata contains the extraction protocol used for each sample, as described in GSE76983.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link[SummarizedExperiment]{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/grun-hsc}.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Grun D et al. (2016). 
#' De novo prediction of stem cell identity using single-cell transcriptome data. 
#' \emph{Cell Stem Cell} 19(2), 266-277. 
#'
#' @examples
#' sce <- GrunHSCData()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData rowData<- colData<-
#' @importFrom S4Vectors DataFrame
GrunHSCData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("grun-bone_marrow-2016", "2023-12-14", realize.assays=TRUE)

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("grun-hsc", version), has.rowdata=FALSE, has.coldata=FALSE)

        # Cleaning up the row data.
        symbol <- sub("__.*", "", rownames(sce))
        loc <- sub(".*__", "", rownames(sce))
        rowData(sce) <- DataFrame(symbol=symbol, chr=loc)

        # Cleaning up the col data.
        cn <- colnames(sce)
        sample <- sub("_.*", "", cn)
        protocol <- ifelse(sample %in% c("JC4", "JC48P2", "JC48P4", "JC48P6", "JC48P7"),
            "sorted hematopoietic stem cells", "micro-dissected cells")
        colData(sce) <- DataFrame(sample=sample, protocol=protocol, row.names=cn)
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rowData(sce)$symbol,
        ensembl=ensembl,
        location=location)
}
