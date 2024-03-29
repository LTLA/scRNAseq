#' Survey of dataset metadata
#'
#' Metadata survey for all available datasets in the \pkg{scRNAseq} package.
#'
#' @param cache,overwrite Arguments to pass to \code{\link{fetchMetadataDatabase}}.
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
#' @seealso
#' \code{\link{searchDatasets}}, to search on the metadata for specific datasets.
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom gypsum cacheDirectory fetchMetadataDatabase
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
surveyDatasets <- function(cache=cacheDirectory(), overwrite=FALSE, latest=TRUE) {
    bpath <- fetchMetadataDatabase(cache=cache, overwrite=overwrite)
    con <- dbConnect(SQLite(), bpath)
    on.exit(dbDisconnect(con))

    stmt <- "SELECT json_extract(metadata, '$') AS meta, versions.asset AS asset, versions.version AS version, path";
    if (!latest) {
        stmt <- paste0(stmt, ", versions.latest AS latest")
    }
    stmt <- paste0(stmt, " FROM paths LEFT JOIN versions ON paths.vid = versions.vid WHERE versions.project = 'scRNAseq'")
    if (latest) {
        stmt <- paste0(stmt, " AND versions.latest = 1")
    }
    everything <- dbGetQuery(con, stmt)

    sanitize_query_to_output(everything, latest)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom jsonlite fromJSON
sanitize_query_to_output <- function(results, latest, meta.name="meta") {
    path <- results$path
    has.slash <- grepl("/", path)
    path[!has.slash] <- NA_character_
    path[has.slash] <- sub("/[^/]+$", "", path[has.slash])
    df <- DataFrame(name = results$asset, version = results$version, path = path)
    if (!latest) {
        df$latest <- results$latest == 1
    }

    all_meta <- lapply(results[[meta.name]], fromJSON, simplifyVector=FALSE)
    df$object <- extract_atomic_from_json(all_meta, function(x) x$applications$takane$type, "character") 
    df$title <- extract_atomic_from_json(all_meta, function(x) x$title, "character")
    df$description <- extract_atomic_from_json(all_meta, function(x) x$title, "character")
    df$taxonomy_id <- extract_charlist_from_json(all_meta, function(x) x$taxonomy_id)
    df$genome <- extract_charlist_from_json(all_meta, function(x) x$genome)

    df$rows <- extract_atomic_from_json(all_meta, function(x) x$applications$takane$summarized_experiment$rows, "integer")
    df$columns <- extract_atomic_from_json(all_meta, function(x) x$applications$takane$summarized_experiment$columns, "integer")
    df$assays <- extract_charlist_from_json(all_meta, function(x) x$applications$takane$summarized_experiment$assays)
    df$column_annotations <- extract_charlist_from_json(all_meta, function(x) x$applications$takane$summarized_experiment$column_annotations)
    df$reduced_dimensions <- extract_charlist_from_json(all_meta, function(x) x$applications$takane$single_cell_experiment$reduced_dimensions)
    df$alternative_experiments <- extract_charlist_from_json(all_meta, function(x) x$applications$takane$single_cell_experiment$alternative_experiments)

    df$bioconductor_version < extract_atomic_from_json(all_meta, function(x) x$bioconductor_version, "character")
    df$maintainer_name < extract_atomic_from_json(all_meta, function(x) x$maintainer_name, "character")
    df$maintainer_email < extract_atomic_from_json(all_meta, function(x) x$maintainer_email, "character")

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
    df$sources <- as(sources, "CompressedList")

    df
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
