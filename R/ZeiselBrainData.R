#' Obtain the Zeisel brain data
#'
#' Obtain the mouse brain single-cell RNA-seq dataset from Zeisel et al. (2015).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Row data contains a single \code{"featureType"} field describing the type of each feature
#' (endogenous genes, mitochondrial genes, spike-in transcripts and repeats).
#' Spike-ins and repeats are stored as separate entries in the \code{\link{altExps}}.
#'
#' Column metadata is provided in the same form as supplied in \url{http://linnarssonlab.org/cortex/}.
#' This contains information such as the cell diameter and the published cell type annotations.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' Spike-in metadata is added using \code{\link{ERCCSpikeInConcentrations}},
#' with molecule counts computed using a volume of 9 nL per cell at a dilution of 1:20000.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zeisel-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Zeisel A et al. (2015). 
#' Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. 
#' \emph{Science} 347(6226), 1138-42.
#'
#' @examples
#' sce <- ZeiselBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps altExp<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics cbind
ZeiselBrainData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("zeisel-brain", version))

    status <- rowData(sce)$featureType
    status[status=="mito"] <- "endogenous"
    sce <- splitAltExps(sce, status, "endogenous")

    # Decorating spike-ins with row metadata. Probably could be less
    # verbose but this is pretty clear.
    spike.exp <- altExp(sce, "ERCC")
    spikedata <- ERCCSpikeInConcentrations(volume = 9, dilution = 20000)
    spikedata <- spikedata[rownames(spike.exp),]

    rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
    rowData(spike.exp)$featureType <- NULL # redundant field; what else would it be!?
    altExp(sce, "ERCC") <- spike.exp

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Mm",
        ensembl=ensembl,
        location=location)
}
