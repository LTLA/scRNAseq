% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CampbellBrainData.R
\name{CampbellBrainData}
\alias{CampbellBrainData}
\title{Obtain the Campbell brain data}
\usage{
CampbellBrainData(ensembl = FALSE, location = TRUE, legacy = FALSE)
}
\arguments{
\item{ensembl}{Logical scalar indicating whether the row names of the returned object should contain Ensembl identifiers.}

\item{location}{Logical scalar indicating whether genomic coordinates should be returned.}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
}
\description{
Obtain the mouse brain single-cell RNA-seq data from Campbell et al. (2017).
}
\details{
Column metadata is provided in the same form as supplied in GSE93374.
This contains information such as the diet of the mice, sex and proposed cell type for each cell.

If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.

If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
Note that this is only performed if \code{ensembl=TRUE}.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/campbell-brain}.
}
\examples{
sce <- CampbellBrainData()

}
\references{
Campbell R et al. (2017).
A molecular census of arcuate hypothalamus and median eminence cell types. 
\emph{Nat. Neurosci.} 20, 484-496.
}
\author{
Aaron Lun
}
