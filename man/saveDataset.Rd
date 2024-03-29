% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveDataset.R
\name{saveDataset}
\alias{saveDataset}
\title{Save a dataset to disk}
\usage{
saveDataset(x, path, metadata)
}
\arguments{
\item{x}{A \linkS4class{SummarizedExperiment} or one of its subclasses.}

\item{path}{String containing the path to a new directory in which to save \code{x}.
Any existing directory is removed before saving \code{x}.}

\item{metadata}{Named list containing metadata for this dataset,
see the schema returned by \code{\link{fetchMetadataSchema}()}.
Note that the \code{applications.takane} property will be automatically added by this function and does not have to be supplied.}
}
\value{
\code{x} and its metadata are saved into \code{path}, and \code{NULL} is invisibly returned.
}
\description{
Save a single-cell dataset to disk, usually in preparation for upload.
}
\examples{
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=matrix(rpois(1000, lambda=1), 100, 10)))
rownames(sce) <- sprintf("GENE_\%i", seq_len(nrow(sce)))
colnames(sce) <- head(LETTERS, 10)

meta <- list(
    title="My dataset",
    description="This is my dataset",
    taxonomy_id="10090",
    genome="GRCh38",
    sources=list(list(provider="GEO", id="GSE12345")),
    maintainer_name="Shizuka Mogami",
    maintainer_email="mogami.shizuka@765pro.com"
)

tmp <- tempfile()
saveDataset(sce, tmp, meta)
list.files(tmp, recursive=TRUE)
alabaster.base::readObject(tmp)

}
\seealso{
\url{https://github.com/ArtifactDB/bioconductor-metadata-index}, on the expected schema for the metadata.

\code{\link{polishDataset}}, to polish \code{x} before saving it.

\code{\link{uploadDirectory}}, to upload the saved contents.
}
\author{
Aaron Lun
}
