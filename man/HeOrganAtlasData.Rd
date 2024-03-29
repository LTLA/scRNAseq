% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeOrganAtlasData.R
\name{HeOrganAtlasData}
\alias{HeOrganAtlasData}
\title{Obtain the He organ atlas data}
\usage{
HeOrganAtlasData(
  tissue = c("Bladder", "Blood", "Common.bile.duct", "Esophagus", "Heart", "Liver",
    "Lymph.node", "Marrow", "Muscle", "Rectum", "Skin", "Small.intestine", "Spleen",
    "Stomach", "Trachea"),
  ensembl = FALSE,
  location = TRUE,
  legacy = FALSE
)
}
\arguments{
\item{tissue}{Character vector specifying the tissues to return.}

\item{ensembl}{Logical scalar indicating whether the output row names should contain Ensembl identifiers.}

\item{location}{Logical scalar indicating whether genomic coordinates should be returned.}

\item{legacy}{Logical scalar indicating whether to pull data from ExperimentHub.
By default, we use data from the gypsum backend.}
}
\value{
A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
}
\description{
Obtain the human cortex single-nuclei RNA-seq data from Hu et al. (2017).
}
\details{
Column data contains the tissue of origin, a variety of per-cell QC metrics well as some cell type annotations.
The reclustered annotations required some assembly:
\itemize{
\item \code{reclustered.broad} was generated based on whether the barcode was present in each \code{*_meta.data.txt} file at \url{https://github.com/bei-lab/scRNA-AHCA}.
\item For each barcode that was present in one of those files, \code{reclustered.fine} was generated based on the label in the \code{annotation} field inside that file.
}

If multiple tissues are requested, counts are only reported for the intersection of genes across all tissues.
This is because the gene annotation in the original count matrices differs across tissues.

If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.

If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
Note that this is only performed if \code{ensembl=TRUE}.

All data are downloaded from ExperimentHub and cached for local re-use.
Specific resources can be retrieved by searching for \code{scRNAseq/he-organ-atlas}.
}
\examples{
if (.Machine$sizeof.pointer > 4) { # too large for 32-bit machines!
    sce <- HeOrganAtlasData()
}

}
\references{
He S et al. (2020).
Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs. 
\emph{Genome Biol} 21, 1:294.
}
\author{
Aaron Lun
}
