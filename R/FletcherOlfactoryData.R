#' Obtain the Fletcher Olfactory data
#'
#' Obtain the mouse olfactory epithelial HBC stem cell differentiation dataset from Fletcher et al. (2017).
#'
#' @param filtered Logical scalar indicating whether to filter out cells that were not used by the authors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details 
#' Column metadata is scraped from GEO, using both the author-supplied \dQuote{phenoData} per-cell annotations and the author-supplied \dQuote{protocolData} per-cell annotations. 
#' The former includes information about the animals and the instruments, while the latter contains QC statistics.
#' 
#' We also included the clustering results from the authors' analysis.
#'
#' If \code{filtered=TRUE}, only the cells used by the authors in their cluster analysis are returned. 
#' Otherwise, the cells not used by the authors will have \code{NA} in the clustering columns of the \code{\link{colData}}.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object. 
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/fletcher-olfactory}.
#'
#' @return 
#' A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Davide Risso
#'
#' @references
#' Fletcher R et al. (2017). 
#' Deconstructing olfactory stem cell trajectories at single-cell resolution. 
#' \emph{Cell Stem Cell} 20, 817-30.
#'
#' @examples
#' sce <- FletcherOlfactoryData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps 
#' @importFrom SummarizedExperiment rowData 
FletcherOlfactoryData <- function(filtered=TRUE, ensembl=FALSE, location=TRUE) {
    version <- "2.6.0"

    sce <- .create_sce(file.path("fletcher-olfactory", version), has.rowdata=TRUE)

    if (filtered) {
        sce <- sce[,sce$retained]
        sce$retained <- NULL
    }

    type <- ifelse(rowData(sce)$Transcript_Type == "brainProjectControl", "control", "endogenous")
    sce <- splitAltExps(sce, type, ref="endogenous")
    
    .convert_to_ensembl(sce, 
        symbols=rowData(sce)$Gene_Symbol, 
        species="Mm",
        ensembl=ensembl,
        location=location)
}
