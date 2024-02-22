#' Obtain the Zilionis lung cancer data
#'
#' Obtain the human/mouse lung cancer single-cell RNA-seq data from Zilionis et al. (2019).
#'
#' @param which String specifying the species to get data for.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param filter Logical scalar indicating if the filtered subset should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata is provided and contains information on the library, donor ID/animal ID, replicate and tissue.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#' 
#' If \code{filter=TRUE}, only cells that have been used in the original analysis are returned.
#' The cells used are specified in the \code{Used} column of the \code{\link{colData}}.
#' 
#' The \code{\link{reducedDim}} contains coordinates of SPRING representations.
#' This may be filled with \code{NA}s for SPRING coordinates computed on a subset of cells (specified in \code{\link{colData}}).
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zilionis-lung}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Jens Preussner
#'
#' @references
#' Zilionis R et al. (2019). 
#' Single-cell transcriptomics of human and mouse lung cancers reveals 
#' conserved myeloid populations across individuals and species. 
#' \emph{Immunity} 50(5), 1317-1334.
#'
#' @examples
#' sce.human <- ZilionisLungData()
#'
#' sce.mouse <- ZilionisLungData("mouse")
#' 
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment reducedDim<- reducedDims<- 
ZilionisLungData <- function(which=c("human", "mouse"), ensembl=FALSE, location=TRUE, filter=FALSE, legacy=FALSE) {
    which <- match.arg(which)

    if (!legacy) {
        sce <- fetchDataset("zilionis-lung-2018")

    } else {
        version <- "2.4.0"
        sce <- .create_sce(file.path("zilionis-lung", version), has.rowdata=FALSE, suffix=which)

        if (filter) {
            sce <- sce[,sce$Used]
        }

        cd <- colData(sce)
        cn <- colnames(cd)
        if (which=="human") {
            which.x <- grep("^x_", cn)
            which.y <- grep("^y_", cn)
            grouping.x <- sub("^x_", "", cn[which.x])
            grouping.y <- sub("^y_", "", cn[which.y])
            stopifnot(identical(grouping.x, grouping.y)) # sanity check.

            all.dims <- list()
            for (i in seq_along(grouping.x)) {
                all.dims[[paste0("SPRING_", grouping.x[i])]] <- cbind(x=cd[,which.x[i]], y=cd[,which.y[i]]) 
            }
            reducedDims(sce) <- all.dims
        } else {
            which.x <- which(cn=="x")
            which.y <- which(cn=="y")
            reducedDim(sce, "SPRING") <- cbind(x=cd[,which.x], y=cd[,which.y])
        } 

        colData(sce) <- colData(sce)[,-c(which.x, which.y)]
    }

    .convert_to_ensembl(sce, 
        ensembl=ensembl,
        species=if (which=="human") "Hs" else "Mm",
        symbols=rownames(sce),
        location=location)
}
