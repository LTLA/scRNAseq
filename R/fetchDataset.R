#' Fetch a dataset from the gypsum backend
#'
#' Fetch a dataset (or its metadata) from the gypsum backend.
#'
#' @param name String containing the name of the dataset.
#' @param version String containing the version of the dataset.
#' @param subpath String containing the path to a subdataset, if \code{name} contains multiple datasets.
#' @param package String containing the name of the package.
#' @param cache,overwrite Arguments to pass to \code{\link[gypsum]{saveVersion}} or \code{\link[gypsum]{saveFile}}.
#' @param ... Further arguments to pass to \code{\link{readObject}}.
#'
#' @return \code{fetchDataset} returns the dataset as a \linkS4class{SummarizedExperiment} or one of its subclasses.
#'
#' \code{fetchMetadata} returns a named list of metadata for the specified dataset.
#'
#' @seealso
#' \code{\link{createMetadata}}, on some of the metadata requirements.
#'
#' \code{\link{saveDataset}} and \code{\link{uploadDirectory}}, to save and upload a dataset.
#'
#' \code{\link{listAvailableVersions}} and friends, to get options for \code{name} and \code{version}.
#' 
#' @author Aaron Lun
#' @examples
#' fetchDataset("zeisel-brain-2015", "2.17.1")
#' fetchMetadata("zeisel-brain-2015", "2.17.1")
#'
#' @export
#' @importFrom alabaster.base altReadObjectFunction readObject
fetchDataset <- function(name, version, subpath=NULL, package="scRNAseq", cache=NULL, overwrite=FALSE, ...) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }
    vpath <- gypsum::saveVersion(package, name, version, cache=cache, overwrite=overwrite)
    provenance <- list(name=name, version=version, package=package, root=vpath)

    opath <- vpath
    if (!is.null(subpath)) {
        subpath <- gsub("/*$", "", subpath)
        opath <- file.path(vpath, subpath)
    }

    old <- altReadObjectFunction(scLoadObject)
    on.exit(altReadObjectFunction(old))
    readObject(opath, scrnaseqmatrix.provenance=provenance, ...)
}

#' @export
fetchMetadata <- function(name, version, path=NULL, package="scRNAseq", cache=NULL, overwrite=FALSE) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }

    if (is.null(path)) {
        path <- "_bioconductor.json"
    } else {
        path <- paste0(path, "/_bioconductor.json")
    }

    mpath <- gypsum::saveFile(package, name, version, path, cache=cache, overwrite=overwrite)
    jsonlite::fromJSON(mpath, simplifyVector=FALSE)
}

#' @importFrom alabaster.base readObjectFile
scLoadObject <- function(path, metadata=NULL, scrnaseqmatrix.provenance=NULL, ...) {
    if (is.null(metadata)) {
        metadata <- readObjectFile(path)
    }
    ans <- readObject(path, metadata=metadata, scrnaseqmatrix.provenance=scrnaseqmatrix.provenance, ...)

    # Heuristic to determine whether something is array-like or not.
    if (is(ans, "DelayedArray") || metadata$type %in% c()) {
        prov <- scrnaseqmatrix.provenance

        # Need to figure out the path inside the version.
        relative.path <- character(0)
        dpath <- path
        while (dirname(dpath) != prov$root) {
            relative.path <- c(basename(dpath), relative.path)
            dpath <- dirname(dpath)
        }

        ans <- ScrnaseqArray(
            name=prov$name, 
            version=prov$version, 
            path=paste(relative.path, collapse="/"), 
            cached=path, 
            package=prov$package, 
            seed=ans
        )
    }

    ans
}
