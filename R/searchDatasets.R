#' Search dataset metadata
#'
#' Search for datasets of interest based on matching text in the associated metadata.
#'
#' @param query String or a \code{gypsum.search.object}, see Examples.
#' @inheritParams surveyDatasets
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
#' searchDatasets("brain")[,c("name", "title")]
#' searchDatasets(defineTextQuery("Neuro%", partial=TRUE))[,c("name", "title")]
#' searchDatasets(defineTextQuery("10090", field="taxonomy_id"))[,c("name", "title")]
#' searchDatasets(
#'    defineTextQuery("GRCm38", field="genome") &
#'    (defineTextQuery("neuro%", partial=TRUE) | 
#'     defineTextQuery("pancrea%", partial=TRUE))
#' )[,c("name", "title")]
#' 
#' @seealso
#' \code{\link{surveyDatasets}}, to easily obtain a listing of all available datasets.
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom gypsum cacheDirectory fetchMetadataDatabase searchMetadataTextFilter
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
searchDatasets <- function(query, cache=cacheDirectory(), overwrite=FALSE, latest=TRUE) {
    filter <- searchMetadataTextFilter(query)

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
    if (!is.null(filter)) {
        stmt <- paste0(stmt, " AND ", filter$where)
        everything <- dbGetQuery(con, stmt, params=filter$parameters)
    } else {
        everything <- dbGetQuery(con, stmt)
    }

    sanitize_query_to_output(everything, latest)
}

#' @export
gypsum::defineTextQuery
