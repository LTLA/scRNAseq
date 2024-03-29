#' Fetch a dataset from the gypsum backend
#'
#' Fetch a dataset (or its metadata) from the gypsum backend.
#'
#' @param name String containing the name of the dataset.
#' @param version String containing the version of the dataset.
#' @param path String containing the path to a subdataset, if \code{name} contains multiple datasets.
#' Defaults to \code{NA} if no subdatasets are present.
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
#' \url{https://github.com/ArtifactDB/bioconductor-metadata-index}, on the expected schema for the metadata.
#'
#' \code{\link{saveDataset}} and \code{\link{uploadDirectory}}, to save and upload a dataset.
#'
#' \code{\link{surveyDatasets}} and \code{\link{listVersions}}, to get possible values for \code{name} and \code{version}.
#' 
#' @author Aaron Lun
#' @examples
#' fetchDataset("zeisel-brain-2015", "2023-12-14")
#' fetchMetadata("zeisel-brain-2015", "2023-12-14")
#'
#' @export
#' @importFrom gypsum cacheDirectory saveVersion
#' @importFrom alabaster.base altReadObjectFunction altReadObject
fetchDataset <- function(name, version, path=NA, package="scRNAseq", cache=cacheDirectory(), overwrite=FALSE, realize.assays=FALSE, realize.reduced.dims=TRUE, ...) {
    version_path <- saveVersion(package, name, version, cache=cache, overwrite=overwrite)

    obj_path <- version_path
    if (!is.na(path)) {
        obj_path <- file.path(version_path, gsub("/*$", "", path))
    }

    old <- altReadObjectFunction(scLoadObject)
    on.exit(altReadObjectFunction(old))
    altReadObject(obj_path, scRNAseq.realize.assays=realize.assays, scRNAseq.realize.reduced.dims=realize.reduced.dims, ...)
}

#' @export
#' @rdname fetchDataset
#' @importFrom jsonlite fromJSON
#' @importFrom gypsum cacheDirectory saveFile
fetchMetadata <- function(name, version, path=NA, package="scRNAseq", cache=cacheDirectory(), overwrite=FALSE) {
    remote_path <- "_bioconductor.json"
    if (!is.na(path)) {
        remote_path <- paste0(path, "/", remote_path)
    }

    local_path <- saveFile(package, name, version, remote_path, cache=cache, overwrite=overwrite)
    fromJSON(local_path, simplifyVector=FALSE)
}

#' @importFrom alabaster.base readObjectFile readObject
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
scLoadObject <- function(path, metadata=NULL, scRNAseq.realize.assays=FALSE, scRNAseq.realize.reduced.dims=TRUE, ...) {
    if (is.null(metadata)) {
        metadata <- readObjectFile(path)
    }
    ans <- readObject(
        path, 
        metadata=metadata, 
        scRNAseq.realize.assays=scRNAseq.realize.assays, 
        scRNAseq.realize.reduced.dims=scRNAseq.realize.reduced.dims, 
        ...
    )

    if (is(ans, "SummarizedExperiment")) {
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

#' @importClassesFrom alabaster.matrix ReloadedArray
#' @importFrom DelayedArray is_sparse type
#' @importClassesFrom Matrix lgCMatrix dgCMatrix
realize_array <- function(x) {
    if (is(x, "ReloadedArray")) {
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
