#' Obtain the La Manno brain data
#'
#' Obtain the mouse/human brain scRNA-seq data from La Manno et al. (2016).
#'
#' @param which A string specifying which dataset should be obtained.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' Column metadata is provided in the same form as supplied in the supplementary tables in GSE71585.
#' This contains information such as the time point and cell type.
#'
#' The various settings of \code{which} will obtain different data sets.
#' \itemize{
#' \item \code{"human-es"}, human embryonic stem cells.
#' \item \code{"human-embryo"}, human embryo midbrain.
#' \item \code{"human-ips"}, human induced pluripotent stem cells.
#' \item \code{"mouse-adult"}, mouse adult dopaminergic neurons.
#' \item \code{"mouse-embryo"}, mouse embryo midbrain.
#' }
#' Unfortunately, each of these datasets uses a different set of features.
#' If multiple datasets are to be used simultaneously, users will have to decide how to merge them,
#' e.g., by taking the intersection of common features across all datasets.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/lamanno-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' La Manno A et al. (2016).
#' Molecular diversity of midbrain development in mouse, human, and stem cells. 
#' \emph{Cell} 167(2), 566-580.
#'
#' @examples
#' sce.h.es <- LaMannoBrainData()
#'
#' sce.h.em <- LaMannoBrainData("human-embryo")
#'
#' sce.h.ip <- LaMannoBrainData("human-ips")
#'
#' sce.m.ad <- LaMannoBrainData("mouse-adult")
#'
#' sce.m.em <- LaMannoBrainData("mouse-embryo")
#' 
#' @export
LaMannoBrainData <- function(which=c("human-es", "human-embryo", "human-ips", "mouse-adult", "mouse-embryo"),
    ensembl=FALSE, location=TRUE) 
{
    version <- "2.0.0"
    sce <- .create_sce(file.path("lamanno-brain", version), has.rowdata=FALSE, suffix=match.arg(which))
    colnames(sce) <- colData(sce)[[grep("Cell_ID", colnames(colData(sce)), ignore.case=TRUE)]]

    .convert_to_ensembl(sce, 
        species=if (grepl("human", which)) "Hs" else "Mm",
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
