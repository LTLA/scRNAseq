% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listVersions.R
\name{listVersions}
\alias{listVersions}
\alias{fetchLatestVersion}
\title{List available versions}
\usage{
listVersions(name)

fetchLatestVersion(name)
}
\arguments{
\item{name}{String containing the name of the dataset.}
}
\value{
For \code{listVersions}, a character vector containing the names of the available versions of the dataset.

For \code{fetchLatestVersion}, a string containing the name of the latest version.
}
\description{
List the available and latest versions for a dataset.
}
\examples{
listVersions("zeisel-brain-2015")
fetchLatestVersion("zeisel-brain-2015")

}
\author{
Aaron Lun
}
