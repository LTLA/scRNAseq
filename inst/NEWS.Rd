\name{scRNAseqnews}
\title{scRNAseq News}
\encoding{UTF-8}

\section{Version 2.4.0}{\itemize{
\item Added the Zilionis Lung dataset (Jens Preussner).

\item Added the Hermann Spermatogenesis dataset (Charlotte Soneson).

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
