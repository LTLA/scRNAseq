#' List available versions
#'
#' List the available and latest versions for a dataset.
#'
#' @param name String containing the name of the dataset.
#'
#' @return 
#' For \code{listVersions}, a character vector containing the names of the available versions of the dataset.
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
