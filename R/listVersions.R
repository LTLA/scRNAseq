#' List available versions
#'
#' List the available and latest versions for an asset (i.e., datasets or collections of datasets).
#'
#' @param name String containing the name of the asset.
#'
#' @return 
#' For \code{listVersions}, a character vector containing the names of the available versions of the \code{name} asset.
#'
#' For \code{fetchLatestVersion}, a string containing the name of the latest version.
#'
#' @author Aaron Lun
#'
#' @examples
#' listVersions("zeisel-brain-2015")
#' fetchLatestVersion("zeisel-brain-2015")
#'
#' @export
listVersions <- function(name) {
    gypsum::listVersions("scRNAseq", name)
}

#' @export
#' @rdname listVersions
#' @importFrom gypsum fetchLatest
fetchLatestVersion <- function(name) {
    fetchLatest("scRNAseq", name)
}
