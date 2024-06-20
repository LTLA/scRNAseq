\name{scRNAseqnews}
\title{scRNAseq News}
\encoding{UTF-8}

\section{Version 2.20.0}{\itemize{
\item Support complex queries in \code{searchDatasets()} via a human-friendly syntax.
}}

\section{Version 2.18.0}{\itemize{
\item Switched to the ArtifactDB representations for the underlying files.
This uses language-agnostic formats (e.g., HDF5, JSON) instead of RDS files to store the various parts of each SingleCellExperiment.
The user experience should be more or less the same as the datasets are indistinguishable once loaded into memory.

\item Added the \code{fetchDataset}, to create SingleCellExperiment objects from the ArtifactDB file representations.
This uses the \pkg{alabaster.base} package to do the loading, with some optional realization of the assays into memory.
Advanced users can achieve faster loading times by keeping the assays as file-backed matrices.

\item Introduced \code{saveDataset} and related functions to facilitate user uploads of their own datasets.
This is accompanied by some step-by-step instructions in the vignette, plus some maintainer instructions in the README.

\item Added \code{searchDatasets} to perform text searches on the metadata for each dataset,
using the SQLite database compiled from the \pkg{gypsum} backend where the files are stored.

\item Updated some datasets to reflect upstream changes (e.g., in ArrayExpress).
Currently, this affects mostly \code{SegerstolpePancreasData}, where ArrayExpress decided to change the names and contents of various column annotations.

\item Soft-deprecation of some redundant pieces of information in each dataset.
Some examples are the column names of the Zilionis data, which were not unique and had no meaning;
or the symbols of the \code{rowData} in the Segerstolpe data, which were redundant with the row names.
These changes will only take effect when \code{fetchDataset} is used directly;
the per-dataset getter functions have appropriate back-compatibility patches to restore this information.

\item All getters now have a \code{legacy=} option to pull RDS files from ExperimentHub instead of the new formats from gypsum.

\item Fixed phase assignments in the Leng ESC dataset to correctly report G2 phase.
}}

\section{Version 2.6.0}{\itemize{
\item Added the Bacher T cell dataset.

\item Added the Bhaduri organoid dataset.

\item Added the Darmanis brain dataset.

\item Added the Ernst spermatogenesis dataset.

\item Added the Fletcher olfactory dataset.

\item Added the Giladi HSC dataset.

\item Added the He organ atlas dataset.

\item Added the Jessa brain dataset.

\item Added the Nowakowski cortex dataset.

\item Added the Pollen glia dataset.

\item Added the Zeisel nervous system dataset.

\item Added the Zhao immune liver dataset.

\item Added the Zhong prefrontal cortex dataset.

\item Added the Bunis HSPC dataset (Dan Bunis).
}}

\section{Version 2.4.0}{\itemize{
\item Added the Zilionis lung dataset (Jens Preussner).

\item Added the Hermann spermatogenesis dataset (Charlotte Soneson).

\item Added the Mair and Kotliarov PBMC datasets (Stephany Orjuela).

\item Added the Stoeckius cell hashing dataset.

\item Added the Wu kidney snRNA-seq dataset.

\item Added the Hu cortex snRNA-seq dataset.

\item Added spike-in concentrations to the \code{altExp} \code{rowData} for various datasets (Alan O'Callaghan).
}}

\section{Version 2.2.0}{\itemize{
\item Removed deprecated inbuilt datasets.

\item Added the Paul HSC dataset.

\item Added a dedicated getter for the ERCC concentrations (Alan O'Callaghan).
}}

\section{Version 2.0.0}{\itemize{
\item Added lots of new ExperimentHub datasets, inspired by simpleSingleCell use cases and Martin Hemberg's website.

\item All outputs are now SingleCellExperiment instances with spike-ins stored as alternative experiments.

\item Deprecated inbuilt datasets in favor of ExperimentHub equivalents.
}}

\section{Version 1.0.0}{\itemize{
\item New package scRNAseq, providing single-cell RNA sequencing datasets.
}}
