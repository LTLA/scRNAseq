#' Obtain the He organ atlas data
#'
#' Obtain the human cortex single-nuclei RNA-seq data from Hu et al. (2017).
#'
#' @param tissue Character vector specifying the tissues to return.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column data contains the tissue of origin, a variety of per-cell QC metrics well as some cell type annotations.
#' The reclustered annotations required some assembly:
#' \itemize{
#' \item \code{reclustered.broad} was generated based on whether the barcode was present in each \code{*_meta.data.txt} file at \url{https://github.com/bei-lab/scRNA-AHCA}.
#' \item For each barcode that was present in one of those files, \code{reclustered.fine} was generated based on the label in the \code{annotation} field inside that file.
#' }
#'
#' If multiple tissues are requested, counts are only reported for the intersection of genes across all tissues.
#' This is because the gene annotation in the original count matrices differs across tissues.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/he-organ-atlas}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' He S et al. (2020).
#' Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs. 
#' \emph{Genome Biol} 21, 1:294.
#'
#' @examples
#' if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
#'     sce <- HeOrganAtlasData()
#' }
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
HeOrganAtlasData <- function(tissue=c("Bladder", "Blood", "Common.bile.duct", "Esophagus", 
    "Heart", "Liver", "Lymph.node", "Marrow", "Muscle", "Rectum", "Skin", "Small.intestine", 
    "Spleen", "Stomach", "Trachea"), ensembl=FALSE, location=TRUE, legacy=FALSE)
{
    if (!legacy) {
        all.sce <- lapply(tissue, function(x) fetchDataset("he-organs-2020", "2023-12-21", path=gsub("\\.", "_", tolower(x)), realize.assays=TRUE))
        if (length(all.sce) > 1) {
            common <- Reduce(intersect, lapply(all.sce, rownames))
            for (i in seq_along(all.sce)) {
                all.sce[[i]] <- all.sce[[i]][common,]
            }
        }
        sce <- do.call(cbind, all.sce)

    } else {
        version <- "2.6.0"
        tissue <- match.arg(tissue, several.ok=TRUE)
        hub <- .ExperimentHub()
        host <- file.path("he-organ-atlas", version)

        all.sce <- lapply(tissue, function(x) .create_sce(host, has.rowdata=FALSE, suffix=x, hub=hub))
        if (length(all.sce) > 1) {
            common <- Reduce(intersect, lapply(all.sce, rownames))
            for (i in seq_along(all.sce)) {
                all.sce[[i]] <- all.sce[[i]][common,]
            }
        }
        sce <- do.call(cbind, all.sce)

        # Pulling down reduced dimensions.
        reddim <- hub[hub$rdatapath==file.path("scRNAseq", host, "reddim.rds")][[1]]
        reddim <- lapply(reddim, function(x) x[colnames(sce),,drop=FALSE])
        reducedDims(sce) <- reddim
    }

    .convert_to_ensembl(sce, 
        species="Hs", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
