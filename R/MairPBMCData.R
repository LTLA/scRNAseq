#' Obtain the Mair CITE-seq data
#'
#' Obtain the Mair PBMC targeted CITE-seq data from Mair et al. (2020).
#'
#' @param mode Character vector specifying whether to return either or both the RNA and ADT counts.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' Column metadata contains the donor identity and cartridge of origin.
#' Some libraries may also be classified as multiplets or have undeterminate origins after hash tag debarcoding.
#'
#' If \code{ensembl=TRUE}, the gene symbols in the RNA data are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link[SummarizedExperiment]{rowRanges}} for the RNA data.
#' Note that this is only performed if \code{ensembl=TRUE}.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/mair-pbmc}.
#' 
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with a single matrix of UMI counts corresponding to the first \code{mode},
#' with an optional alternative Experiment if there is a second \code{mode}.
#'
#' @author
#' Stephany Orjuela,
#' with modifications from Aaron Lun
#'
#' @references
#' Mair C et al. (2020).
#' A targeted multi-omic analysis approach measures protein expression and low-abundance transcripts on the single-cell level.
#' \emph{Cell Rep.} 31, 107499
#'
#' @examples
#' sce <- MairPBMCData()
#' 
#' @export
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SummarizedExperiment colData<- rowData
#' @importFrom SingleCellExperiment SingleCellExperiment altExps 
MairPBMCData <- function(mode=c("rna", "adt"), ensembl=FALSE, location=TRUE, legacy=FALSE) {
    mode <- match.arg(mode, c("rna", "adt"), several.ok=TRUE)

    if (!legacy && identical(mode, c("rna", "adt"))) {
        sce <- fetchDataset("mair-pbmc-2020", "2024-04-18", realize.assays=TRUE)

    } else {
        version <- "2.4.0"
        tag <- "mair-pbmc"
        hub <- ExperimentHub()

        collated <- list()
        for (x in mode) {
            collated[[x]] <- .create_sce(file.path(tag, version), hub=hub, 
                has.rowdata=TRUE, has.coldata=FALSE, suffix=x)
        }

        if ("rna" %in% names(collated)) {
            collated[["rna"]] <- .convert_to_ensembl(collated[["rna"]],
                symbols=rowData(collated[["rna"]])$Symbol,
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
