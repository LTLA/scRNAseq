#' Fetch a dataset from the gypsum backend
#'
#' Fetch a dataset (or its metadata) from the gypsum backend.
#'
#' @param name String containing the name of the dataset.
#' @param version String containing the version of the dataset.
#' @param path String containing the path to a subdataset, if \code{name} contains multiple datasets.
#' Defaults to \code{"."} if no subdatasets are present.
#' @param package String containing the name of the package.
#' @param cache,overwrite Arguments to pass to \code{\link[gypsum]{saveVersion}} or \code{\link[gypsum]{saveFile}}.
#' @param realize.assays,realize.reduced.dims Logical scalars indicating whether to realize assays and reduced dimensions into memory.
#' Dense and sparse \linkS4class{ReloadedArray} objects are converted into ordinary arrays and \linkS4class{dgCMatrix} objects, respectively.
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
#' fetchDataset("zeisel-brain-2015", "2023-12-14")
#' fetchMetadata("zeisel-brain-2015", "2023-12-14")
#'
#' @export
#' @importFrom alabaster.base altReadObjectFunction altReadObject
fetchDataset <- function(name, version, path=".", package="scRNAseq", cache=NULL, overwrite=FALSE, realize.assays=FALSE, realize.reduced.dims=TRUE, ...) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }
    version_path <- gypsum::saveVersion(package, name, version, cache=cache, overwrite=overwrite)
    provenance <- list(name=name, version=version, package=package, root=normalizePath(version_path))

    obj_path <- version_path
    if (path != ".") {
        obj_path <- file.path(version_path, gsub("/*$", "", path))
    }

    old <- altReadObjectFunction(scLoadObject)
    on.exit(altReadObjectFunction(old))
    altReadObject(obj_path, scRNAseq.array.provenance=provenance, scRNAseq.realize.assays=realize.assays, scRNAseq.realize.reduced.dims=realize.reduced.dims, ...)
}

#' @export
fetchMetadata <- function(name, version, path=".", package="scRNAseq", cache=NULL, overwrite=FALSE) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }

    if (path == ".") {
        remote_path <- "_bioconductor.json"
    } else {
        remote_path <- paste0(path, "/_bioconductor.json")
    }

    local_path <- gypsum::saveFile(package, name, version, remote_path, cache=cache, overwrite=overwrite)
    jsonlite::fromJSON(local_path, simplifyVector=FALSE)
}

#' @importFrom alabaster.base readObjectFile readObject
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
scLoadObject <- function(path, metadata=NULL, scRNAseq.array.provenance=NULL, scRNAseq.realize.assays=FALSE, scRNAseq.realize.reduced.dims=TRUE, ...) {
    if (is.null(metadata)) {
        metadata <- readObjectFile(path)
    }
    ans <- readObject(
        path, 
        metadata=metadata, 
        scRNAseq.array.provenance=scRNAseq.array.provenance, 
        scRNAseq.realize.assays=scRNAseq.realize.assays, 
        scRNAseq.realize.reduced.dims=scRNAseq.realize.reduced.dims, 
        ...
    )

    if (is(ans, "ReloadedArray")) {
        prov <- scRNAseq.array.provenance 

        # Need to figure out the relative path inside the version. This assumes
        # that 'prov$root' was also normalized so that paths are comparable.
        relative.path <- character(0)
        dpath <- normalizePath(path)
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
            seed=ans@seed@seed
        )

    } else if (is(ans, "SummarizedExperiment")) {
        if (scRNAseq.realize.assays) {
            for (y in assayNames(ans)) {
                assay(ans, y, withDimnames=FALSE) <- realize_array(assay(ans, y, withDimnames=FALSE))
            }
        }

        if (is(ans, "SingleCellExperiment")) {
            if (scRNAseq.realize.reduced.dims) {
                for (z in reducedDimNames(ans)) {
                    reducedDim(ans, z, withDimnames=FALSE) <- realize_array(reducedDim(ans, z, withDimnames=FALSE))
                }
            }
        }
    }

    ans
}

#' @importFrom DelayedArray is_sparse type
realize_array <- function(x) {
    if (is(x, "DelayedArray")) {
        if (is_sparse(x)) {
            if (type(x) == "logical") {
                x <- as(x, "lgCMatrix")
            } else {
                x <- as(x, "dgCMatrix")
            }
        } else {
            x <- as.array(x)
        }
    }
    x
}
