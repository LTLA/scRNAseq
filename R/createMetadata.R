#' Helper function to create dataset metadata
#'
#' Create metadata for a dataset, typically to be used in \code{\link{saveDataset}}.
#' This can also be created directly via the \code{\link{list}} constructor, for users who already know all the required fields.
#'
#' @param title String containing the title for this dataset.
#' @param description String containing the description for this dataset.
#' @param taxonomy Character vector of NCBI taxonomy IDs for the species involved in this dataset.
#' @param genome Character vector of Ensembl genome builds used in this dataset.
#' For example, this should be \code{"GRCm38"} instead of the UCSC notation \code{"mm10"}.
#' Some genomes may not yet be added to the schema - in such cases, contact the package maintainer.
#' @param sources A list of lists, where each inner list specifies a contributing source of the dataset.
#' Each inner list should be named with the \code{provider} and \code{id} elements.
#' \itemize{
#' \item For \code{provider="GEO"}, the \code{id} should be formatted as \dQuote{GSExxx}.
#' \item For \code{provider="ArrayExpress"}, the \code{id} should be formatted as \dQuote{E-MTAB-xxx}.
#' \item For \code{provider="PubMed"}, the \code{id} should only contain digits.
#' \item For \code{provider="DOI"}, the \code{id} should be a DOI without the \dQuote{doi:} prefix.
#' \item For \code{provider="other"}, the \code{id} should be a URL.
#' }
#' @param maintainer.name String containing the name of the dataset maintainer, i.e., the person who prepared \code{x}.
#' This should contain space-separated first, middle and last names.
#' @param maintainer.email String containing the email for the maintainer.
#' @param version String containing the Bioconductor version.
#' This defaults to the version in the current R session.
#'
#' @return A named list of metadata.
#'
#' @examples
#' meta <- createMetadata(
#'     title="My dataset",
#'     description="This is my dataset",
#'     taxonomy.id=c("10090", "9606"),
#'     genome=c("GRCh38", "GRCm39"),
#'     sources=list(
#'         list(provider="GEO", id="GSE12345"),
#'         list(provider="PubMed", id="67890")
#'     ),
#'     maintainer.name="Shizuka Mogami",
#'     maintainer.email="mogami.shizuka@765pro.com"
#' )
#'
#' # Just makes a regular list:
#' str(meta)
#'
#' # Comparing it against the schema:
#' gypsum::validateMetadata(meta)
#'
#' @author Aaron Lun
#'
#' @export
createMetadata <- function(title, description, taxonomy.id, genome, sources, maintainer.name, maintainer.email, bioc.version=NULL) {
    if (is.null(bioc.version)) {
        bioc.version <- as.character(BiocManager::version())
    }
    list(
        title=title,
        description=description,
        taxonomy_id=taxonomy.id,
        genome=genome,
        sources=sources,
        maintainer_name=maintainer.name,
        maintainer_email=maintainer.email,
        bioconductor_version=bioc.version
    )
}
