% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NowakowskiCortexData.R
\name{NowakowskiCortexData}
\alias{NowakowskiCortexData}
\title{Obtain the Nowakowski cortex data}
\usage{
NowakowskiCortexData(ensembl = FALSE, location = TRUE, legacy = FALSE)
}
\arguments{
\item{ensembl}{Logical scalar indicating whether the output row names should contain Ensembl identifiers.}

\item{location}{Logical scalar indicating whether genomic coordinates should be returned.}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of TPMs.
The \code{\link{reducedDims}} contains an assortment of dimensionality reduction results.
}
\description{
Obtain the human cortex single-cell RNA-seq dataset from Nowakowski et al. (2017).
}
\details{
Column metadata includes the presumed cell type (\code{WGCNAcluster}), patient and tissue region of origin. 
A variety of dimensionality reduction results are also provided.

If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.

If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
This is only performed when \code{ensembl=TRUE}.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/nowakowski-cortex}.
}
\examples{
sce <- NowakowskiCortexData()

}
\references{
Nowakowski S et al. (2017). 
Spatiotemporal gene expression trajectories reveal developmental hierarchies of the human cortex. 
\emph{Science} 358, 1318-1323.
}
\author{
Aaron Lun
}
