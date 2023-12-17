#' List available datasets
#'
#' Summary information for all available datasets in the \pkg{scRNAseq} package.
#'
#' @param cache,overwrite Arguments to pass to \code{\link[gypsum]{fetchMetadataDatabase}}.
#'
#' @return 
#' A \linkS4class{DataFrame} where each row corresponds to a dataset, containing various pieces of metadata.
#' Some columns may be lists to capture 1:many mappings.
#'
#' @author Aaron Lun
#'
#' @examples
#' listAvailableDatasets()
#' 
#' @export
#' @importFrom S4Vectors DataFrame
listAvailableDatasets <- function(cache=NULL, overwrite=FALSE) {
    if (is.null(cache)) {
        cache <- gypsum::cacheDirectory()
    }
    bpath <- gypsum::fetchMetadataDatabase("bioconductor", cache=cache, overwrite=overwrite)
    spath <- gypsum::fetchMetadataDatabase("scRNAseq", cache=cache, overwrite=overwrite)

    bcon <- DBI::dbConnect(RSQLite::SQLite(), bpath)
    bcore <- DBI::dbReadTable(bcon, "core", check.names=FALSE)
    bfree <- DBI::dbReadTable(bcon, "free_text", check.names=FALSE)
    btaxid <- DBI::dbReadTable(bcon, "multi_taxonomy_id", check.names=FALSE)
    bgenome <- DBI::dbReadTable(bcon, "multi_genome", check.names=FALSE)
    bsources <- DBI::dbReadTable(bcon, "multi_sources", check.names=FALSE)

    scon <- DBI::dbConnect(RSQLite::SQLite(), spath)
    score <- DBI::dbReadTable(scon, "core", check.names=FALSE)
    sass <- DBI::dbReadTable(scon, "multi_assay_names", check.names=FALSE)
    sred <- DBI::dbReadTable(scon, "multi_reduced_dimension_names", check.names=FALSE)
    salt <- DBI::dbReadTable(scon, "multi_alternative_experiment_names", check.names=FALSE)

    keys <- bcore$`_key`
    m_free2core <- match(keys, bfree$`_key`)
    m_s2core <- match(keys, bfree$`_key`)

    DataFrame(
        asset = bcore$`_asset`,
        version = bcore$`_version`,
        path = bcore$`_path`,
        object = bcore$`_object`,

        title = bfree$title[m_free2core],
        description = bfree$description[m_free2core],
        taxonomy_id = splitAsList(btaxid$item, factor(btaxid$`_key`, keys)),
        genome = splitAsList(bgenome$item, factor(bgenome$`_key`, keys)),
        sources = I(splitAsList(DataFrame(bsources[,c("provider", "id")]), factor(bsources$`_key`, keys))),

        nrow = score$nrow[m_s2core],
        ncol = score$ncol[m_s2core],
        assay_names = splitAsList(sass$item, factor(sass$`_key`, keys)),
        reduced_dimension_names = splitAsList(sred$item, factor(sred$`_key`, keys)),
        alternative_experiment_names = splitAsList(salt$item, factor(salt$`_key`, keys)),

        bioconductor_version = bcore$bioconductor_version,
        maintainer_name = bfree$maintainer_name[m_free2core],
        maintainer_email = bcore$maintainer_email
    )
}
