#' Reprocessed single-cell data sets
#'
#' Obtain the legacy count matrices for three publicly available single-cell RNA-seq datasets.
#' Raw sequencing data were downloaded from NCBI's SRA or from EBI's ArrayExpress,
#' aligned to the relevant genome build and used to quantify gene expression.
#'
#' @param assays Character vector specifying one or more assays to return.
#' Choices are \code{"tophat_counts"}, \code{"cufflinks_fpkm"}, \code{"rsem_counts"} and \code{"rsem_tpm"}.
#' If \code{NULL}, all assays are returned.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @return
#' A \linkS4class{SingleCellExperiment} object containing one or more expression matrices of counts and/or TPMs,
#' depending on \code{assays}.
#' 
#' @details
#' \code{ReprocessedFluidigmData} returns a dataset of 65 human neural cells from Pollen et al. (2014), 
#' each sequenced at high and low coverage (SRA accession SRP041736).
#' 
#' \code{ReprocessedTh2Data} returns a dataset of 96 mouse T helper cells from Mahata et al. (2014),
#' obtained from ArrayExpress accession E-MTAB-2512.
#' Spike-in counts are stored in the \code{"ERCC"} entry of the \code{\link{altExps}}.
#'
#' \code{ReprocessedAllenData} return a dataset of 379 mouse brain cells from Tasic et al. (2016).
#' This is a re-processed subset of the data from \code{\link{TasicBrainData}},
#' and contains spike-in information stored as in the \code{\link{altExps}}.
#'
#' In each dataset, the first columns of the \code{colData} are sample quality metrics from FastQC and Picard.
#' The remaining fields were obtained from the original study in their GEO/SRA submission
#' and/or as Supplementary files in the associated publication.
#' These two categories of \code{colData} are distinguished by a \code{which_qc} element in the \code{\link{metadata}},
#' which contains the names of the quality-related columns in each object.
#' 
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/legacy-allen},
#' \code{scRNAseq/legacy-fluidigm} or \code{scRNAseq/legacy-th2}.
#'
#' @section Pre-processing details:
#' FASTQ files were either obtained directly from ArrayExpress, 
#' or converted from SRA files (downloaded from the Sequence Read Archive) using the SRA Toolkit.
#'
#' Reads were aligned with TopHat (v. 2.0.11) to the appropriate reference genome (GRCh38 for human samples, GRCm38 for mouse). 
#' RefSeq mouse gene annotation (GCF_000001635.23_GRCm38.p3) was downloaded from NCBI on Dec. 28, 2014. 
#' RefSeq human gene annotation (GCF_000001405.28) was downloaded from NCBI on Jun. 22, 2015.
#'
#' featureCounts (v. 1.4.6-p3) was used to compute gene-level read counts.
#' Cufflinks (v. 2.2.0) was used to compute gene-leve FPKMs.
#' Reads were also mapped to the transcriptome using RSEM (v. 1.2.19) to compute read counts and TPM's.
#' 
#' FastQC (v. 0.10.1) and Picard (v. 1.128) were used to compute sample quality control (QC) metrics. 
#' However, no filtering on the QC metrics has been performed for any dataset.
#'
#' @references
#' Pollen AA et al. (2014). 
#' Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. 
#' \emph{Nat. Biotechnol.} 32(10), 1053-8.
#'
#' Mahata B et al. (2014).
#' Single-cell RNA sequencing reveals T helper cells synthesizing steroids de novo to contribute to immune homeostasis. 
#' \emph{Cell Rep}, 7(4), 1130-42.
#'
#' Tasic A et al. (2016). 
#' Adult mouse cortical cell taxonomy revealed by single cell transcriptomics.
#' \emph{Nat. Neurosci.} 19(2), 335-46.
#'
#' @examples
#' sce <- ReprocessedAllenData()
#'
#' @export
#' @rdname ReprocessedData
#' @importFrom SingleCellExperiment splitAltExps
#' @importFrom SummarizedExperiment assays assays<-
ReprocessedAllenData <- function(assays=NULL, ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("legacy", "2023-12-18", "allen", realize.assays=TRUE)
        if (!is.null(assays)) {
            assays(sce, withDimnames=FALSE) <- assays(sce, withDimnames=FALSE)[assays]
        }
    } else {
        version <- "1.10.0"
        sce <- .create_sce_legacy(file.path("legacy-allen", version), assays)
        status <- ifelse(grepl("^ERCC-[0-9]+$", rownames(sce)), "ERCC", "endogenous")
        sce <- splitAltExps(sce, status, ref="endogenous")
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}

#' @export
#' @rdname ReprocessedData
#' @importFrom SingleCellExperiment splitAltExps
ReprocessedTh2Data <- function(assays=NULL, ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("legacy", "2023-12-18", "th2", realize.assays=TRUE)
        if (!is.null(assays)) {
            assays(sce, withDimnames=FALSE) <- assays(sce, withDimnames=FALSE)[assays]
        }
    } else {
        version <- "1.10.0"
        sce <- .create_sce_legacy(file.path("legacy-th2", version), assays)
        status <- ifelse(grepl("^ERCC-[0-9]+$", rownames(sce)), "ERCC", "endogenous")
        sce <- splitAltExps(sce, status, ref="endogenous")
    }

    .convert_to_ensembl(sce, 
        species="Mm", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}

#' @export
#' @rdname ReprocessedData
ReprocessedFluidigmData <- function(assays=NULL, ensembl=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("legacy", "2023-12-18", "fluidigm", realize.assays=TRUE)
        if (!is.null(assays)) {
            assays(sce, withDimnames=FALSE) <- assays(sce, withDimnames=FALSE)[assays]
        }
    } else {
        version <- "1.10.0"
        sce <- .create_sce_legacy(file.path("legacy-fluidigm", version), assays)
    }

    .convert_to_ensembl(sce, 
        species="Hs", 
        symbols=rownames(sce),
        ensembl=ensembl,
        location=location)
}
