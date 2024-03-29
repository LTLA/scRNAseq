% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BachMammaryData.R
\name{BachMammaryData}
\alias{BachMammaryData}
\title{Obtain the Bach mammary data}
\usage{
BachMammaryData(
  samples = c("NP_1", "NP_2", "G_1", "G_2", "L_1", "L_2", "PI_1", "PI_2"),
  location = TRUE,
  legacy = FALSE
)
}
\arguments{
\item{samples}{A character vector with at least one element, specifying which samples(s) to retrieve.}

\item{location}{Logical scalar indicating whether genomic coordinates should be returned.}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
}
\description{
Obtain the mouse mammary gland single-cell RNA-seq data from Bach et al. (2017).
}
\details{
Column metadata is extracted from the sample annotation in GSE106273,
and refers to the developmental stage of the mammary gland.

If multiple samples are specified in \code{samples}, the count matrices will be \code{cbind}ed together.
Cells originating from different samples are identifiable by the \code{"Sample"} field in the column metadata.

If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/bach-mammary}.
}
\examples{
sce <- BachMammaryData(samples="NP_1")

}
\references{
Bach K et al. (2017).
Differentiation dynamics of mammary epithelial cells revealed by single-cell RNA sequencing. 
\emph{Nat Commun.} 8(1), 2128
}
\author{
Aaron Lun
}
