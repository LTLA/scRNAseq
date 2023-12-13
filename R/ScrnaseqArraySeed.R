#' scRNAseq-derived array 
#'
#' An array that was loaded as part of a dataset from \code{\link{fetchDataset}}.
#' This allows methods to avoid re-saving the array and uploading a duplicate copy when, e.g., updating an existing dataset.
#'
#' @param name String containing the name of the dataset containing the array.
#' Alternatively an existing ScrnaseqArraySeed, which is returned without modification.
#' @param version String containing the version of the dataset containing the array.
#' @param package String containing the name of the package containing the dataset.
#' @param path String containing the path to the array inside the dataset.
#' @param cached String containing the path to the locally cached representation of the array.
#' @param seed Contents of the loaded array.
#' If \code{NULL}, this is created by calling \code{\link{readObject}} on \code{cached}.
#' @param ... Further arguments to pass to \code{\link{readObject}} when \code{seed=NULL}.
#'
#' @return
#' For the constructors, an instance of the \linkS4class{ScrnaseqArraySeed} or \linkS4class{ScrnaseqArray}.
#'
#' @details
#' The ScrnaseqArraySeed is a subclass of the \linkS4class{WrapperArraySeed} and will just forward all operations to the underlying \code{seed}.
#' Its main purpose is to track the gypsum project-asset-version combination that was originally used to generate \code{seed}, 
#' which enables optimizations during \code{\link{saveDataset}} and \code{\link{uploadDirectory}}.
#'
#' @name ScrnaseqArraySeed
#' @aliases
#' ScrnaseqArraySeed-class
#' ScrnaseqArray-class
#' ScrnaseqMatrix-class
#' DelayedArray,ScrnaseqArraySeed-method
#' matrixClass,ScrnaseqArray-method
#' coerce,ScrnaseqArray,ScrnaseqMatrix-method
#' coerce,ScrnaseqMatrix,ScrnaseqArray-method
#' path,ScrnaseqArraySeed-method
#' saveObject,ScrnaseqArray-method
#'
#' @export
ScrnaseqArraySeed <- function(name, version, path, cached, package="scRNAseq", seed=NULL, ...) {
    if (is(name, "ScrnaseqArraySeed")) {
        return(path)
    }
    if (is.null(seed)) {
        seed <- readObject(cached, ...)
    }
    while (is(seed, "DelayedArray")) {
        seed <- seed@seed
    }
    new("ScrnaseqArraySeed", name=name, version=version, path=path, cached=cached, package=package, seed=seed)
}

#' @export
#' @rdname ScrnaseqArraySeed
#' @importFrom DelayedArray DelayedArray
ScrnaseqArray <- function(path, seed=NULL, ...) {
    DelayedArray(ScrnaseqArraySeed(path, seed=seed, ...))
}

#' @export
#' @importClassesFrom alabaster.matrix WrapperArraySeed
setClass("ScrnaseqArraySeed", contains="WrapperArraySeed", slots=c(name="character", version="character", path="character", cached="character", package="character"))

#' @export
#' @importClassesFrom DelayedArray DelayedArray
setClass("ScrnaseqArray", contains="DelayedArray", slots=c(seed = "ScrnaseqArraySeed"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ScrnaseqMatrix", contains=c("ScrnaseqArray", "DelayedMatrix"))

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ScrnaseqArraySeed", function(seed) new_DelayedArray(seed, Class="ScrnaseqArray"))

#' @export
#' @importFrom DelayedArray matrixClass
setMethod("matrixClass", "ScrnaseqArray", function(x) "ScrnaseqMatrix")

# Overrides copied from DelayedArray::ConstantArray.
#' @importFrom S4Vectors new2
setAs("ScrnaseqArray", "ScrnaseqMatrix", function(from) new2("ScrnaseqMatrix", from))
setAs("ScrnaseqMatrix", "ScrnaseqArray", function(from) from)

#' @export
#' @importFrom alabaster.base saveObject
setMethod("saveObject", "ScrnaseqArray", function(x, path, reloadedarray.reuse.files="link", ...) {
    s <- x@seed

    manifest <- list.files(s@cached, recursive=TRUE)
    dir.create(path)

    for (y in manifest) {
        combined <- file.path(path, y)
        dir.create(dirname(combined), showWarnings=FALSE, recursive=TRUE)
        if (!file.symlink(file.path(s@cached, y), file.path(path, y))) {
            stop("failed to create symbolic link '", y, "' from '", s@cached, "' to '", path, "'")
        }
    }

    write(file=".link", toJSON(list(project=s@package, asset=s@name, version=s@version, path=s@path), auto_unbox=TRUE))
    invisible(NULL)
})
