#' Search dataset metadata
#'
#' Search for datasets of interest based on matching text in the associated metadata.
#'
#' @param query String containing a query in a human-readable syntax or a \link{gypsum.search.clause}, see Examples.
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
#' searchDatasets("Neuro%")[,c("name", "title")]
#' searchDatasets("taxonomy_id:10090")[,c("name", "title")]
#' searchDatasets("(genome: GRCm38 AND neuro%) OR pancrea%")[,c("name", "title")]
#'
#' # We can also use gypsum search clauses via the gsc() function:
#' searchDatasets(gsc(asset="he-organs-2020"))[,c("path")]
#' searchDatasets(
#'     (gsc(asset="%brain%", partial=TRUE) |
#'     gsc(asset="%neur%", partial=TRUE) |
#'     gsc(text="%neur%", partial=TRUE) |
#'     gsc(text="%brain%", partial=TRUE)) &
#'     gsc(field="taxonomy_id", text="10090")
#' )[,c("name", "title")]
#' 
#' @seealso
#' \code{\link{surveyDatasets}}, to easily obtain a listing of all available datasets.
#'
#' \code{\link{translateTextQuery}}, for details on the human-readable query syntax.
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom gypsum cacheDirectory fetchMetadataDatabase searchMetadataFilter translateTextQuery
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
searchDatasets <- function(query, cache=cacheDirectory(), overwrite=FALSE, latest=TRUE) {
    if (is.character(query)) {
        query <- translateTextQuery(query)
    }
    filter <- searchMetadataFilter(query)

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
#' @importFrom gypsum defineTextQuery
gypsum::defineTextQuery

#' @export
#' @importFrom gypsum gsc
gypsum::gsc
