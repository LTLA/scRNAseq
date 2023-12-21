#' Save a dataset to disk
#'
#' Save a single-cell dataset to disk, usually in preparation for upload.
#'
#' @param x A \linkS4class{SummarizedExperiment} or one of its subclasses.
#' @param path String containing the path to a new directory in which to save \code{x}.
#' Any existing directory is removed before saving \code{x}.
#' @param metadata Named list containing metadata, usually generated using \code{\link{createMetadata}} or \code{\link{fetchMetadata}}.
#'
#' @return \code{x} and its metadata are saved into \code{path}, and \code{NULL} is invisibly returned.
#'
#' @seealso
#' \code{\link{createMetadata}}, to create the metadata.
#'
#' \code{\link{polishDataset}}, to polish \code{x} before saving it.
#'
#' \code{\link{uploadDirectory}}, to upload the saved contents.
#' 
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(list(counts=matrix(rpois(1000, lambda=1), 100, 10)))
#' rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
#' colnames(sce) <- head(LETTERS, 10)
#'
#' meta <- createMetadata(
#'     title="My dataset",
#'     description="This is my dataset",
#'     taxonomy.id="10090",
#'     genome="GRCh38",
#'     sources=list(list(provider="GEO", id="GSE12345")),
#'     maintainer.name="Shizuka Mogami",
#'     maintainer.email="mogami.shizuka@765pro.com"
#' )
#' 
#' tmp <- tempfile()
#' saveDataset(sce, tmp, meta)
#' list.files(tmp, recursive=TRUE)
#' alabaster.base::readObject(tmp)
#' 
#' @export
#' @importFrom alabaster.base saveObject
#' @importMethodsFrom alabaster.sce saveObject
saveDataset <- function(x, path, metadata) {
    metadata$taxonomy_id <- I(metadata$taxonomy_id)
    metadata$genome <- I(metadata$genome)
    contents <- jsonlite::toJSON(metadata, pretty=4, auto_unbox=TRUE)
    gypsum::validateMetadata(contents, schema=gypsum::fetchMetadataSchema())

    unlink(path, recursive=TRUE)
    alabaster.base::saveObject(x, path)

    write(contents, file=file.path(path, "_bioconductor.json"))
    invisible(NULL)
}
