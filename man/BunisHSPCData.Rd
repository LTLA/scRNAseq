% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BunisHSPCData.R
\name{BunisHSPCData}
\alias{BunisHSPCData}
\title{Obtain the Bunis haematopoietic stem and progenitor cell data}
\usage{
BunisHSPCData(filtered = TRUE, legacy = FALSE)
}
\arguments{
\item{filtered}{Logical scalar or "cells" indicating whether to:
\itemize{
\item \code{TRUE}: filter out cells that were not used by the authors.
\item \code{"cells"}: filter out empty droplets as filtered out by cell ranger.
\item \code{FALSE}: no filtering
}}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
}
\description{
Obtain the human fetal, newborn, and adult haematopoietic stem and progenitor cell single-cell RNA-seq dataset from Bunis et al. (2021).
}
\details{
Column metadata is recreated from GEO using the author-supplied TSV of per-cell annotations, or retrieved from a processed version of the data shared by authors via figshare.
This contains information such as the tissue & sample of origin, age group, likely cell type, and Developmental Stage Scoring.
Within DevStageScoring element of the column metadata are the applied results ('<cell_type>_scores') of random forest regression trained on the fetal (score = 0) and adult (score = 1) cells of individual cell types indicated by ('<cell_type>_inTraining').

If \code{filtered=TRUE}, only the cells used by the authors in their final analysis are returned.
Otherwise, an additional \code{retained} field will be present in the \code{\link{colData}}, indicating whether the cell was retained by the authors.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/bunis-hspc}.
}
\examples{
sce <- BunisHSPCData()

}
\references{
Bunis DG et al. (2021). 
Single-Cell Mapping of Progressive Fetal-to-Adult Transition in Human Naive T Cells
\emph{Cell Rep.} 34(1): 108573
}
\author{
Daniel Bunis
}
