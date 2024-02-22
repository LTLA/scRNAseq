#' Obtain the Zhong prefrontal cortex data
#'
#' Obtain the human prefrontal cortex single-cell RNA-seq dataset from Zhong et al. (2018).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' Column metadata is scraped from GEO and includes week of gestation, gender and likely cell type.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' This is only performed when \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zhong-prefrontal}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Zhong S et al. (2018). 
#' A single-cell RNA-seq survey of the developmental landscape of the human prefrontal cortex. 
#' \emph{Nature} 555, 524-528.
#'
#' @examples
#' sce <- ZhongPrefrontalData()
#' 
#' @export
ZhongPrefrontalData <- function(ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("zhong-prefrontal-2018", "2023-12-22", realize.assays=TRUE)
    } else {
        version <- "2.6.0"
        sce <- .create_sce(file.path("zhong-prefrontal", version), has.rowdata=FALSE)
    }

    .convert_to_ensembl(sce, 
        symbols=rownames(sce), 
        species="Hs",
        ensembl=ensembl,
        location=location)
}
