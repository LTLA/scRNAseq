#' Obtain the Segerstolpe pancreas data
#'
#' Download the human pancreas single-cell RNA-seq (scRNA-seq) dataset from Segerstolpe et al. (2016)
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Row data contains fields for the gene symbol and RefSeq transcript IDs corresponding to each gene.
#' The rows of the output object are named with the symbol, but note that these are not unique.
#'
#' Column metadata were extracted from the \code{Characteristics} fields of the SDRF file for ArrayExpress E-MTAB-5061.
#' This contains information such as the cell type labels and patient status.
#'
#' Count data for ERCC spike-ins are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#' Estimated numbers of spike-in molecules are provided in the \code{\link{rowData}} of this entry. 
#' Note that these concentrations are incorrect for donor H1, as 100 uL of spike-in mixture were added for this donor, rather than 25 uL for all others.
#' 
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
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
SegerstolpePancreasData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("segerstolpe-pancreas-2016", "2023-12-19", realize.assays=TRUE)

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("segerstolpe-pancreas", version))
        rownames(sce) <- rowData(sce)$symbol

        status <- ifelse(grepl("^ERCC-[0-9]+", rowData(sce)$refseq), "ERCC", "endogenous")
        sce <- splitAltExps(sce, status, ref="endogenous")
        
        ## This is wrong for one donor - donor "H1" has 100ul rather than 25ul
        spike.exp <- altExp(sce, "ERCC")
        spikedata <- ERCCSpikeInConcentrations(volume = 25, dilution = 40000)
        spikedata <- spikedata[rownames(spike.exp), ]
        rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
        altExp(sce, "ERCC") <- spike.exp
    }

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
