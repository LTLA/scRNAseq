% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchDatasets.R
\name{searchDatasets}
\alias{searchDatasets}
\title{Search dataset metadata}
\usage{
searchDatasets(
  query,
  cache = cacheDirectory(),
  overwrite = FALSE,
  latest = TRUE
)
}
\arguments{
\item{query}{String containing a query in a human-readable syntax or a \link[gypsum]{gypsum.search.clause}, see Examples.}

\item{cache, overwrite}{Arguments to pass to \code{\link[gypsum]{fetchMetadataDatabase}}.}

\item{latest}{Whether to only consider the latest version of each dataset.}
}
\value{
A \link[S4Vectors]{DataFrame} where each row corresponds to a (sub)dataset, containing various columns of metadata.
Some columns may be lists to capture 1:many mappings.
}
\description{
Search for datasets of interest based on matching text in the associated metadata.
}
\details{
The returned DataFrame contains the usual suspects like the title and description for each dataset,
the number of rows and columns, the organisms and genome builds involved,
whether the dataset has any pre-computed reduced dimensions, and so on.
More details can be found in the Bioconductor metadata schema at \url{https://github.com/ArtifactDB/bioconductor-metadata-index}. 

If a dataset contains multiple subdatasets, each subdataset is reported as a separate row in the DataFrame.
The location of subdataset is provided in the \code{path} column.
If a dataset does not contain any subdatasets, the \code{path} entry will be set to \code{NA}.
}
\examples{
searchDatasets("brain")[,c("name", "title")]
searchDatasets("Neuro\%")[,c("name", "title")]
searchDatasets("taxonomy_id:10090")[,c("name", "title")]
searchDatasets("(genome: GRCm38 AND neuro\%) OR pancrea\%")[,c("name", "title")]

# We can also use gypsum search clauses via the gsc() function.
# Here, 'asset' is analogous to the 'name' of the dataset.
searchDatasets(gsc(asset="he-organs-2020"))[,c("path")]
searchDatasets(
    (gsc(asset="\%brain\%", partial=TRUE) |
    gsc(asset="\%neur\%", partial=TRUE) |
    gsc(text="\%neur\%", partial=TRUE) |
    gsc(text="\%brain\%", partial=TRUE)) &
    gsc(field="taxonomy_id", text="10090")
)[,c("name", "title")]

}
\seealso{
\code{\link{surveyDatasets}}, to easily obtain a listing of all available datasets.

\code{\link[gypsum]{translateTextQuery}}, for details on the human-readable query syntax.
}
\author{
Aaron Lun
}
