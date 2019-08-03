#' Obtain the Lun spike-in data
#'
#' Download and cache the Lun spike-in single-cell RNA-seq (scRNA-seq) dataset from ExperimentHub,
#' returning a \linkS4class{SingleCellExperiment} object for further use.
#'
#' @param which String specifying whether the 416B or trophoblast data should be obtained.
#' 
#' @details
#' This function provides the spike-in scRNA-seq data from Lun et al. (2017)
#' in the form of a \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' Row data contains a single \code{"Length"} field describing the total exonic length of each feature.
#'
#' Column metadata is provided in the same form as supplied in E-MTAB-5522.
#' This contains information such as the cell type, plate of origin, spike-in addition order and oncogene induction. 
#'
#' Two sets of spike-ins are available for each dataset - SIRVs and ERCCs - and stored as alternative experiments.
#' 
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @references
#' Lun et al. (2017). 
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
#' @importFrom SingleCellExperiment splitSCEByAlt
LunSpikeInData <- function(which=c("416b", "tropho")) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("lun-spikein", version), suffix=match.arg(which))

    spike.type <- rep("endogenous", nrow(sce))
    spike.type[grep("ERCC", rownames(sce))] <- "ERCC"
    spike.type[grep("SIRV", rownames(sce))] <- "SIRV"

    splitSCEByAlt(sce, spike.type, ref="endogenous")
}
