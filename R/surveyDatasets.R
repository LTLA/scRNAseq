#' Survey of dataset metadata
#'
#' Metadata survey for all available datasets in the \pkg{scRNAseq} package.
#'
#' @param cache,overwrite Arguments to pass to \code{\link[gypsum]{fetchMetadataDatabase}}.
#' @param latest Whether to only consider the latest version of each dataset.
#'
#' @return 
#' A \linkS4class{DataFrame} where each row corresponds to a dataset, containing various columns of metadata.
#' Some columns may be lists to capture 1:many mappings.
#'
#' @details
#' The returned DataFrame contains the usual suspects like the title and description for each dataset,
#' the number of rows and columns, the organisms and genome builds involved,
#' whether the dataset has any pre-computed reduced dimensions, and so on.
#' More details can be found in the Bioconductor metadata schema at \url{https://github.com/ArtifactDB/bioconductor-metadata-index}. 
#'
#' @author Aaron Lun
#'
#' @examples
#' surveyDatasets()
#' 
#' @export
#' @importFrom S4Vectors DataFrame
surveyDatasets <- function(cache=NULL, overwrite=FALSE, latest=TRUE) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }
    bpath <- gypsum::fetchMetadataDatabase(cache=cache, overwrite=overwrite)
    con <- DBI::dbConnect(RSQLite::SQLite(), bpath)
    on.exit(DBI::dbDisconnect(con))

    stmt <- "SELECT json_extract(metadata, '$') AS meta, versions.project AS project, versions.asset AS asset, versions.version AS version, path";
    if (!latest) {
        stmt <- paste0(stmt, ", versions.latest AS latest")
    }
    stmt <- paste0(stmt, " FROM paths LEFT JOIN versions ON paths.vid = versions.vid WHERE versions.project = 'scRNAseq'")
    if (latest) {
        stmt <- paste0(stmt, " AND versions.latest = 1")
    }
    everything <- DBI::dbGetQuery(con, stmt)

    output <- DataFrame(
        asset = everything$asset,
        version = everything$version,
        path = everything$path
    )
    if (!latest) {
        output$latest <- everything$latest == 1
    }

    all_meta <- lapply(everything$meta, jsonlite::fromJSON, simplifyVector=FALSE)
    output$object <- extract_atomic_from_json(all_meta, function(x) x$takane$type, "character") 
    output$title <- extract_atomic_from_json(all_meta, function(x) x$title, "character")
    output$description <- extract_atomic_from_json(all_meta, function(x) x$title, "character")
    output$taxonomy_id <- extract_charlist_from_json(all_meta, function(x) x$taxonomy_id)
    output$genome <- extract_charlist_from_json(all_meta, function(x) x$genome)

    output$rows <- extract_atomic_from_json(all_meta, function(x) x$takane$summarized_experiment$rows, "integer")
    output$columns <- extract_atomic_from_json(all_meta, function(x) x$takane$summarized_experiment$columns, "integer")
    output$assays <- extract_charlist_from_json(all_meta, function(x) x$takane$summarized_experiment$assays)
    output$column_annotations <- extract_charlist_from_json(all_meta, function(x) x$takane$summarized_experiment$column_annotations)
    output$reduced_dimensions <- extract_charlist_from_json(all_meta, function(x) x$takane$single_cell_experiment$reduced_dimensions)
    output$alternative_experiments <- extract_charlist_from_json(all_meta, function(x) x$takane$single_cell_experiment$alternative_experiments)

    output$bioconductor_version < extract_atomic_from_json(all_meta, function(x) x$bioconductor_version, "character")
    output$maintainer_name < extract_atomic_from_json(all_meta, function(x) x$maintainer_name, "character")
    output$maintainer_email < extract_atomic_from_json(all_meta, function(x) x$maintainer_email, "character")

    sources <- vector("list", length(all_meta))
    for (i in seq_along(all_meta)) {
        cursources <- all_meta[[i]]$sources
        if (is.null(cursources)) {
            sources[[i]] <- DataFrame(provider=character(0), id=character(0), version=character(0))
        } else {
            sources[[i]] <- DataFrame(
                provider = extract_atomic_from_json(cursources, function(x) x$provider, "character"),
                id = extract_atomic_from_json(cursources, function(x) x$id, "character"),
                version = extract_atomic_from_json(cursources, function(x) x$version, "character")
            )
        }
    }
    output$sources <- as(sources, "CompressedList")

    output
}

extract_atomic_from_json <- function(metadata, extract, type) {
    vapply(metadata, function(y) {
        x <- extract(y)
        if (is.null(x)) {
            as(NA, type)
        } else {
            x
        }
    }, vector(type, 1))
}

extract_charlist_from_json <- function(metadata, extract) {
    output <- lapply(metadata, function(y) {
        x <- extract(y)
        if (is.null(y)) {
            character(0)
        } else {
            x
        }
    })
    as(output, "CompressedList")
}
