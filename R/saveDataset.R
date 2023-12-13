#' Save a dataset to disk
#'
#' Save a single-cell dataset to disk, usually in preparation for upload.
#'
#' @param x A \linkS4class{SummarizedExperiment} or one of its subclasses.
#' @param path String containing the path to a new directory in which to save \code{x}.
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
#' @return \code{x} and its metadata are saved into \code{path}, and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(list(counts=matrix(rpois(1000, lambda=1), 100, 10)))
#' rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
#' colnames(sce) <- head(LETTERS, 10)
#'
#' tmp <- tempfile()
#' saveDataset(sce, tmp,
#'     title="My dataset",
#'     description="This is my dataset",
#'     taxonomy="10090",
#'     genome="GRCh38",
#'     sources=list(list(provider="GEO", id="GSE12345")),
#'     maintainer.name="Shizuka Mogami",
#'     maintainer.email="mogami.shizuka@765pro.com"
#' )
#'
#' list.files(tmp, recursive=TRUE)
#' alabaster.base::readObject(tmp)
#' 
#' @export
saveDataset <- function(x, path, title, description, taxonomy, genome, sources, maintainer.name, maintainer.email, version=NULL) {
    if (is.null(version)) {
        version <- BiocManager::version()
    }

    metadata <- list(
        title=title,
        description=description,
        taxonomy_id=I(taxonomy),
        genome=I(genome),
        sources=sources,
        maintainer_name=maintainer.name,
        maintainer_email=maintainer.email,
        bioconductor_version=as.character(version)
    )
    contents <- jsonlite::toJSON(metadata, pretty=4, auto_unbox=TRUE)
    gypsum::validateMetadata(contents, schema=gypsum::fetchMetadataSchema())

    alabaster.base::saveObject(x, path)
    write(contents, file=file.path(path, "_bioconductor.json"))
    invisible(NULL)
}
