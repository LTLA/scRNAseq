#' Obtain the Hu cortex data
#'
#' Obtain the mouse cortex single-nuclei RNA-seq data from Hu et al. (2017).
#'
#' @param mode Character vector indicating whether to return data for the 3T3 cells or the mouse cortex.
#' @param samples Character vector indicating whether to return data for specific samples, see Details.
#' If specified, this overrides \code{mode}.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata includes the mode and sample corresponding to each cell/nuclei.
#' Available samples are:
#' \itemize{
#' \item \code{"cell-3T3"} and \code{"nuclei-3T3"}, generated from the 3T3 cell line.
#' \item \code{"nuclei-ctx-X"}, nuclei generated from the cortex of animal number X (from 1 to 13).
#' \item \code{"nuclei-ctx-salineX"} or \code{"nuclei-ctx-PTZX"}, nuclei generated from the cortex of saline- or PTZ-treated mice.
#' X represents the replicate number and can be 1 or 2.
#' }
#'
#' If multiple modes are requested, counts are only reported for the intersection of genes across all modes.
#' This is because the gene annotation in the original count matrices differs across modes.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/wu-kidney}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Hu P et al. (2017). 
#' Dissecting cell-type composition and activity-dependent transcriptional state in mammalian brains by massively parallel single-nucleus RNA-seq.
#' \emph{Mol. cell} 68, 1006-1015.
#'
#' @examples
#' sce <- HuCortexData("3T3")
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
HuCortexData <- function(mode=c("ctx", "3T3"), samples=NULL, ensembl=FALSE, location=TRUE, legacy=FALSE) {
    mode <- match.arg(mode)

    if (!legacy) {
        if (mode == "ctx") {
            mode <- "cortex"
        }
        sce <- fetchDataset("hu-cortex-2017", "2023-12-20", path=mode, realize.assays=TRUE)

    } else {
        version <- "2.4.0"

        avail.samples <- c(
            "cell-3T3", "nuclei-3T3",
            sprintf("nuclei-ctx-%i", seq_len(13)),
            "nuclei-ctx-saline1", "nuclei-ctx-PTZ1", "nuclei-ctx-saline2", "nuclei-ctx-PTZ2"
        )
        avail.modes <- sub("^[^-]+-", "", avail.samples)
        avail.modes <- sub("-.*", "", avail.modes)
        
        if (!is.null(samples)) {
            samples <- match.arg(samples, avail.samples, several.ok=TRUE)
        } else {
            mode <- match.arg(mode, several.ok=TRUE)
            samples <- character(0)
            for (m in mode) { # convoluted to respect ordering of 'mode'.
                samples <- c(samples, avail.samples[avail.modes==m])
            }
        }
        sample.idx <- match(samples, avail.samples)

        all.sce <- list()
        for (s in sample.idx) {
            current <- .create_sce(file.path("hu-cortex", version), has.rowdata=FALSE, has.coldata=FALSE, suffix=avail.samples[s])
            current$Mode <- avail.modes[s]
            current$Sample <- avail.samples[s]
            all.sce <- c(all.sce, list(current))
        }

        if (length(all.sce) > 1) {
            common <- Reduce(intersect, lapply(all.sce, rownames))
            for (i in seq_along(all.sce)) {
                all.sce[[i]] <- all.sce[[i]][common,]
            }
        }

        sce <- do.call(cbind, all.sce)
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
