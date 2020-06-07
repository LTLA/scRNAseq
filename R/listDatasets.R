#' List all available datasets
#'
#' Summary information for all available datasets in the \pkg{scRNAseq} package.
#'
#' @details
#' A study may contribute multiple datasets if they cannot be reasonably combined (e.g., different species).
#' If the authors provided more cells than were used in the analysis, 
#' the total number of cells prior to author filtering is reported.
#'
#' @return 
#' A \linkS4class{DataFrame} where each row corresponds to a dataset, containing the fields:
#' \itemize{
#' \item \code{Reference}, a Markdown-formatted citation to \code{scripts/ref.bib} in the \pkg{scRNAseq} installation directory.
#' \item \code{Taxonomy}, an identifier for the organism.
#' \item \code{Part}, the part of the organism being studied.
#' \item \code{Number}, the total number of cells.
#' \item \code{Call}, the relevant R call required to generate the dataset.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' listDatasets()
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.csv
listDatasets <- function() {
    path <- system.file("extdata", "manifest.csv", package="scRNAseq")
    DataFrame(read.csv(path, stringsAsFactors=FALSE))
}

