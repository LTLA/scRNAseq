% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GrunPancreasData.R
\name{GrunPancreasData}
\alias{GrunPancreasData}
\title{Obtain the Grun pancreas data}
\usage{
GrunPancreasData(ensembl = FALSE, location = TRUE, legacy = FALSE)
}
\arguments{
\item{ensembl}{Logical scalar indicating whether the output row names should contain Ensembl identifiers.}

\item{location}{Logical scalar indicating whether genomic coordinates should be returned.}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
}
\description{
Obtain the human pancreas single-cell RNA-seq data from Grun et al. (2016).
}
\details{
Row metadata contains fields for the symbol and chromosomal location of each gene,
as derived from the row names.

Column metadata is derived from the column names of the count matrix with the sample annotations in GSE81076.
This includes the donor identity for each cell and the type of sample.

The \code{"ERCC"} entry in the \code{\link{altExps}} contains count data for the ERCC spike-in transcripts.

If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.

If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
Note that this is only performed if \code{ensembl=TRUE}.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/grun-pancreas}.
}
\examples{
sce <- GrunPancreasData()

}
\references{
Grun D et al. (2016). 
De novo prediction of stem cell identity using single-cell transcriptome data. 
\emph{Cell Stem Cell} 19(2), 266-277.
}
\author{
Aaron Lun,
using additional metadata obtained by Vladimir Kiselev.
}
