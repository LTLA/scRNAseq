#' Obtain the Wu kidney data
#'
#' Obtain the mouse kidney single-nuclei RNA-seq data from Wu et al. (2019).
#'
#' @param mode String indicating whether to return data for healthy and/or diseased donors.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Column metadata includes the single-cell technology and whether they came from a diseased or healthy individual.
#' 
#' If \code{mode} specifies both healthy and disease donors,
#' counts are only reported for the intersection of genes that are present for both donors.
#' This is because the original count matrices had differences in their annotation.
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
#' Wu H et al. (2019). 
#' Advantages of single-nucleus over single-cell RNA sequencing of adult kidney: rare cell types and novel cell states revealed in fibrosis.
#' \emph{J. Am. Soc. Nephrol.} 30, 23-32.
#'
#' @examples
#' sce <- WuKidneyData("disease")
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
WuKidneyData <- function(mode=c("healthy", "disease"), ensembl=FALSE, location=TRUE) {
    version <- "2.4.0"
    mode <- match.arg(mode, several.ok=TRUE)

    all.sce <- list()
    for (m in mode) {
        current <- .create_sce(file.path("wu-kidney", version), has.rowdata=FALSE, has.coldata=FALSE, suffix=m)

        if (m=="healthy") {
            current$Technology <- sub("_.*", "", colnames(current))
        } else {
            current$Technology <- "sNuc-10x"
        }
        current$Status <- m

        all.sce[[m]] <- current
    }

    if (length(all.sce) > 1) {
        common <- Reduce(intersect, lapply(all.sce, rownames))
        for (i in seq_along(all.sce)) {
            all.sce[[i]] <- all.sce[[i]][common,]
        }
    }

    sce <- do.call(cbind, all.sce)

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
