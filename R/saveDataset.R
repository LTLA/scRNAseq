#' Save a dataset to disk
#'
#' Save a single-cell dataset to disk, usually in preparation for upload.
#'
#' @param x A \linkS4class{SummarizedExperiment} or one of its subclasses.
#' @param path String containing the path to a new directory in which to save \code{x}.
#' Any existing directory is removed before saving \code{x}.
#' @param metadata Named list containing metadata for this dataset,
#' see the schema returned by \code{\link[gypsum]{fetchMetadataSchema}()}.
#'
#' @return \code{x} and its metadata are saved into \code{path}, and \code{NULL} is invisibly returned.
#'
#' @seealso
#' \url{https://github.com/ArtifactDB/bioconductor-metadata-index}, on the expected schema for the metadata.
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
#' meta <- list(
#'     title="My dataset",
#'     description="This is my dataset",
#'     taxonomy_id="10090",
#'     genome="GRCh38",
#'     sources=list(list(provider="GEO", id="GSE12345")),
#'     maintainer_name="Shizuka Mogami",
#'     maintainer_email="mogami.shizuka@765pro.com"
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
    schema <- gypsum::fetchMetadataSchema()

    if (is.null(metadata$bioconductor_version)) {
        metadata$bioconductor_version <- as.character(BiocManager::version())
    }
    metadata$taxonomy_id <- I(metadata$taxonomy_id)
    metadata$genome <- I(metadata$genome)
    gypsum::validateMetadata(metadata, schema) # First validation for user-supplied content.

    unlink(path, recursive=TRUE)
    alabaster.base::saveObject(x, path, ReloadedArray.reuse.files="symlink")

    takane <- list(type = readObjectFile(path)$type)

    if (is(x, "SummarizedExperiment")) {
        takane$summarized_experiment = list(
           rows = nrow(x),
           columns = ncol(x),
           assays = I(assayNames(x)),
           column_annotations = I(colnames(colData(x)))
        )

        if (is(x, "SingleCellExperiment")) {
            takane$single_cell_experiment <- list(
                reduced_dimensions = I(reducedDimNames(x)),
                alternative_experiments = I(altExpNames(x))
            )
        }
    } else if (is(x, "DataFrame")) {
        takane$data_frame <- list(
            rows=nrow(x),
            column_names=I(colnames(x))
        )
    }
    metadata$applications <- c(metadata$applications, list(takane=takane))

    # Second validation with the takane metadata.
    contents <- jsonlite::toJSON(metadata, pretty=4, auto_unbox=TRUE)
    gypsum::validateMetadata(contents, schema=schema)
    write(contents, file=file.path(path, "_bioconductor.json"))
    invisible(NULL)
}
