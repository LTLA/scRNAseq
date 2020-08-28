#' Obtain the Kotliarov CITE-seq data
#'
#' Obtain the Kotliarov PBMC CITE-seq data from Kotliarov et al. (2020).
#'
#' @param mode Character vector specifying whether to return either or both the RNA and ADT counts.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' 
#' 20 samples classified as high and low responders were divided into 2 experimental batches.
#' Each batch contains 5 high and 5 low responders. 
#' The 10 samples per batch are mixed and distributed across the 6 lanes in the plate.
#' 
#' Relevant column metadata contains the following columns: 
#' * sample: Sample identifier with _d0.
#' * sample_id: Sample identifier.
#' * adjmfc.time: type of responder.
#' * tenx_lane: for each batch (H1B) there are 6 lanes (e.g. H1B1ln1, ln2, ln3...).
#' * batch: one of 1:2.
#' * barcode_check: barcode identifier.
#' * hash_ and hto_ columns: HTOdemux outputs.
#' * DEMUXLET. columns: demuxlet outputs.
#' * joint_classification_global: HTOdemux and demuxlet joint classification.
#' * nGene: number of genes as defined in Seurat's CreateSeuratObject().
#' * nUMI: number of UMIs, same as above.
#' * pctMT: percent of mitocondrial reads. 
#' Counts have not been filtered using the above 3 columns.
#'  
#'
#' If \code{ensembl=TRUE}, the gene symbols in the RNA data are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} for the RNA data.
#' Note that this is only performed if \code{ensembl=TRUE}.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/kotliarov-pbmc}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts corresponding to the first \code{mode},
#' with an optional alternative Experiment if there is a second \code{mode}.
#'
#' @author
#' Stephany Orjuela,
#' with modifications from Aaron Lun
#'
#' @references
#' Kotliarov et al. (2020).
#' Broad immune activation underlies shared set point signatures for vaccine responsiveness in healthy individuals and disease activity in patients with lupus.
#' \emph{Nat. Med.} 26, 618–629
#'
#' @examples
#' sce <- KotliarovPBMCData()
#' 
#' @export
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment SingleCellExperiment altExps 
KotliarovPBMCData <- function(mode=c("rna", "adt"), ensembl=FALSE, location=TRUE) {
    mode <- match.arg(mode, c("rna", "adt"), several.ok=TRUE)
    version <- "2.4.0"
    tag <- "kotliarov-pbmc"
    hub <- ExperimentHub()
  
    collated <- list()
    for (x in mode) {
      collated[[x]] <- .create_sce(file.path(tag, version), hub=hub, 
                                   has.rowdata=TRUE, has.coldata=FALSE, suffix=x)
    }
  
    if ("rna" %in% names(collated)) {
      collated[["rna"]] <- .convert_to_ensembl(collated[["rna"]],
                                               symbols=rowData(sce)$SYMBOL,
                                               species="Hs",
                                               ensembl=ensembl,
                                               location=location)
    }
  
    sce <- collated[[1]]
    altExps(sce) <- collated[-1]
    colData(sce) <- hub[hub$rdatapath==file.path("scRNAseq", tag, "coldata.rds")][[1]] 
  
    sce
}
