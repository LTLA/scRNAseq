#' List available versions
#'
#' List the available and latest versions for an asset (i.e., datasets or collections of datasets).
#'
#' @param package String containing the name of the package.
#' @param name String containing the name of the asset.
#'
#' @return 
#' For \code{listAvailableVersion}, a character vector containing the names of the available versions of the \code{name} asset.
#'
#' For \code{fetchLatestVersion}, a string containing the name of the latest version.
#'
#' @author Aaron Lun
#'
#' @examples
#' listAvailableVersions("zeisel-brain-2015")
#' fetchLatestVersion("zeisel-brain-2015")
#'
#' @export
#' @rdname listAvailableAssets
listAvailableVersions <- function(name, package="scRNAseq") {
    gypsum::listVersions(package, name)
}

#' @export
#' @rdname listAvailableAssets
fetchLatestVersion <- function(name, package="scRNAseq") {
    gypsum::fetchLatest(package, name)
}
