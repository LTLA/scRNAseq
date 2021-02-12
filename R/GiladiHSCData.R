#' Obtain the Giladi HSC data
#'
#' Obtain the mouse haematopoietic stem cell single-cell RNA-seq and CRISPR-seq dataset from Giladi et al. (2018).
#'
#' @param filtered Logical scalar indicating whether to filter out cells that were not used by the authors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers, when \code{mode} contains \code{"rna"}.
#' @param location Logical scalar indicating whether genomic coordinates should be returned, when \code{mode} contains \code{"rna"}.
#' @param mode Character vector specifying which modalities should be returned.
#' 
#' @details
#' Column metadata is scraped from GEO using the author-supplied TSV of per-cell annotations. 
#' This contains information such as the batch of origin for each cell plus an array of FACS measurements per cell.
#'
#' If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
#' Otherwise, an additional \code{filtered} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors. 
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#' For row names with multiple semi-colon-delimited symbols, the last symbol is used for matching against the Ensembl annotation.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' This is only relevant when \code{ensembl=TRUE}.
#'
#' If \code{mode} contains multiple modalities, the intersection of cells that are present in both modalities is returned.
#' This is because not all cells have data across both modalities.
#' If \code{mode} contains only one modality, all cells for that modality are returned.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/giladi-hsc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a matrix of UMI counts for the scRNA-seq or CRISPR-seq data.
#' Alternatively, an object with both count matrices, where the second modality is stored as an alternative Experiment.
#'
#' @author Aaron Lun
#'
#' @references
#' Giladi A et al. (2018).
#' Single-cell characterization of haematopoietic progenitors and their trajectories in homeostasis and perturbed haematopoiesis. 
#' \emph{Nat Cell Biol} 20, 836-846
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- GiladiHSCData()
#' }
#' @export
#' @importFrom SingleCellExperiment splitAltExps altExp altExp<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics cbind
GiladiHSCData <- function(mode=c("rna", "crispr"), filtered=TRUE, ensembl=FALSE, location=TRUE) {
    mode <- match.arg(mode, several.ok=TRUE)

    version <- "2.6.0"

    collated <- list()
    if ("rna" %in% mode) {
        sce <- .create_sce(file.path("giladi-hsc", version), has.rowdata=FALSE, suffix="rna")

        if (filtered) {
            sce <- sce[,sce$retained]
            sce$filtered <- NULL
        }

        sce <- .convert_to_ensembl(sce, 
            symbols=sub(".*;", "", rownames(sce)),
            species="Mm",
            ensembl=ensembl,
            location=location)

        collated$rna <- sce
    } 

    if ("crispr" %in% mode) {
        sce <- .create_sce(file.path("giladi-hsc", version), suffix="crispr")
        collated$crispr <- sce
    }

    if (length(collated) > 1) {
        keep <- Reduce(intersect, lapply(collated, colnames))
        for (i in seq_along(collated)) {
            collated[[i]] <- collated[[i]][,keep]
        }
    }

    primary <- collated[[1]]
    altExps(primary) <- collated[-1]
    primary
}
