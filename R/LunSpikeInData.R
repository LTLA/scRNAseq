#' Obtain the Lun spike-in data
#'
#' Obtain the spike-in single-cell RNA-seq data from Lun et al. (2017).
#'
#' @param which String specifying whether the 416B or trophoblast data should be obtained.
#' 
#' @details
#'
#' Row data contains a single \code{"Length"} field describing the total exonic length of each feature.
#'
#' Column metadata is provided in the same form as supplied in E-MTAB-5522.
#' This contains information such as the cell type, plate of origin, spike-in addition order and oncogene induction. 
#'
#' Two sets of spike-ins are available for each dataset - SIRVs and ERCCs - and stored as alternative experiments.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/lun-spikein}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Lun ATL et al. (2017). 
#' Assessing the reliability of spike-in normalization for analyses of single-cell RNA sequencing data. 
#' \emph{Genome Res.} 27(11), 1795-1806.
#'
#' @examples
#' sce <- LunSpikeInData()
#'
#' sce <- LunSpikeInData("tropho")
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment splitAltExps
LunSpikeInData <- function(which=c("416b", "tropho")) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("lun-spikein", version), suffix=match.arg(which))

    spike.type <- rep("endogenous", nrow(sce))
    spike.type[grep("ERCC", rownames(sce))] <- "ERCC"
    spike.type[grep("SIRV", rownames(sce))] <- "SIRV"

    splitAltExps(sce, spike.type, ref="endogenous")
}
