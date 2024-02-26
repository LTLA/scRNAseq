#' Obtain the Kotliarov CITE-seq data
#'
#' Obtain the Kotliarov PBMC CITE-seq data from Kotliarov et al. (2020).
#'
#' @param mode Character vector specifying whether to return either or both the RNA and ADT counts.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' This dataset contains 20 samples from 2 experimental batches, where each batch contains 5 high and 5 low responders. 
#' The 10 samples per batch were mixed and distributed across the 6 lanes using a cell hashing approach.
#' 
#' The column metadata contains the following fields:
#' \itemize{
#' \item \code{sample*}: identifiers for the sample of origin for each cell.
#' \item \code{adjmfc.time}: type of responder for each sample.
#' \item \code{tenx_lane}: 10X lane from which each cell was collected.
#' \item \code{batch}: the batch of origin.
#' \item \code{barcode_check}: barcode identifier.
#' \item \code{hash_*} and \code{hto_*} columns: \pkg{HTOdemux} outputs.
#' \item \code{DEMUXLET.*} columns: \pkg{demuxlet} outputs.
#' \item \code{joint_classification_global}: \pkg{HTOdemux} and \pkg{demuxlet} joint classification.
#' \item \code{nGene}: number of genes as defined from \pkg{Seurat}'s \code{CreateSeuratObject}.
#' \item \code{nUMI}: number of UMIs as defined from \pkg{Seurat}'s \code{CreateSeuratObject}.
#' \item \code{pctMT}: percent of mitochondrial reads as defined from \pkg{Seurat}'s \code{CreateSeuratObject}.
#' }
#' Note, no filtering has been performed based on the quality control metrics.
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
KotliarovPBMCData <- function(mode=c("rna", "adt"), ensembl=FALSE, location=TRUE, legacy=FALSE) {
    mode <- match.arg(mode, c("rna", "adt"), several.ok=TRUE)

    if (!legacy && identical(mode, c("rna", "adt"))) {
        sce <- fetchDataset("kotliarov-pbmc-2020", "2023-12-20", realize.assays=TRUE)

    } else {
        version <- "2.4.0"
        tag <- "kotliarov-pbmc"
        hub <- .ExperimentHub()
        
        collated <- list()
        for (x in mode) {
            collated[[x]] <- .create_sce(file.path(tag, version), hub=hub, 
                has.rowdata=FALSE, has.coldata=FALSE, suffix=x)
        }

        if ("rna" %in% names(collated)) {
            collated[["rna"]] <- .convert_to_ensembl(collated[["rna"]],
                symbols=rownames(collated[["rna"]]),
                species="Hs",
                ensembl=ensembl,
                location=location)
        }
        
        sce <- collated[[1]]
        altExps(sce) <- collated[-1]
        colData(sce) <- hub[hub$rdatapath==file.path("scRNAseq", tag, version, "coldata.rds")][[1]] 
    }

    sce
}
