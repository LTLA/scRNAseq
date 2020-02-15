#' Obtain the Romanov brain data
#'
#' Obtain the mouse brain single-cell RNA-seq dataset from Romanov et al. (2017).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Column metadata is provided in the same form as supplied in GSE74672.
#' This contains information such as the reporter gene expressed in each cell, the mouse line, dissection type and so on.
#'
#' Counts for ERCC spike-ins are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#' Note that some of the spike-in rows have \code{NA} observations for some (but not all) cells.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/romanov-brain}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun,
#' based on code by Vladimir Kiselev and Tallulah Andrews.
#'
#' @references
#' Romanov RA et al. (2017).
#' Molecular interrogation of hypothalamic organization reveals distinct dopamine neuronal subtypes. 
#' \emph{Nat. Neurosci.} 20, 176-188.
#'
#' @examples
#' sce <- RomanovBrainData()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
RomanovBrainData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("romanov-brain", version), has.rowdata=FALSE)

    status <- ifelse(grepl("^ERCC-[0-9]+", rownames(sce)), "ERCC", "endogenous")
    sce <- splitAltExps(sce, status, ref="endogenous")

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
