#' List all paths for a dataset
#'
#' List the available paths to subdatasets within a version of a given dataset.
#' 
#' @inheritParams fetchDataset
#' @param include.metadata Logical scalar indicating whether to report the metadata for each subdataset.
#'
#' @return If \code{include.metadata=FALSE}, a character vector containing the paths to subdatasets within the specified version of the dataset.
#' If \code{name} does not contain any subdatasets, \code{NA} is returned.
#'
#' Otherwise, a \link[S4Vectors]{DFrame} is returned containing the metadata for each subdataset, e.g., the title and description.
#' More details can be found in the Bioconductor metadata schema at \url{https://github.com/ArtifactDB/bioconductor-metadata-index}. 
#'
#' @author Aaron Lun
#'
#' @examples
#' listPaths("he-organs-2020", "2023-12-21")
#' listPaths("he-organs-2020", "2023-12-21", include.metadata=TRUE)
#' listPaths("zeisel-brain-2015", "2023-12-14") # no subdatasets
#' listPaths("zeisel-brain-2015", "2023-12-14", include.metadata=TRUE) 
#'
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom gypsum saveFile listFiles
listPaths <- function(name, version, package="scRNAseq", cache=cacheDirectory(), overwrite=FALSE, include.metadata=FALSE) {
    all.files <- listFiles(package, name, version, include..=FALSE)
    all.files <- all.files[basename(all.files) == "_bioconductor.json"]
    if (!include.metadata) {
        if (identical(all.files, "_bioconductor.json")) {
            return(NA_character_)
        } else {
            return(dirname(all.files))
        }
    }

    metadata <- list()
    for (x in seq_along(all.files)) {
        local.path <- saveFile(package, name, version, path=all.files[x], cache=cache, overwrite=overwrite)
        metadata[[x]] <- paste(readLines(local.path), collapse="\n")
    }

    sanitize_query_to_output(list(
        meta=metadata,
        path=all.files,
        asset=rep(name, length(metadata)),
        version=rep(version, length(metadata))
    ), latest=TRUE)
}
