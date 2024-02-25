#' Obtain the Lun spike-in data
#'
#' Obtain the spike-in single-cell RNA-seq data from Lun et al. (2017).
#'
#' @param which String specifying whether the 416B or trophoblast data should be obtained.
#' @param split.oncogene Logical scalar indicating whether the oncogene should be split to a separate \code{\link{altExp}}.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#' 
#' @details
#' Row data contains a single \code{"Length"} field describing the total exonic length of each feature.
#'
#' Column metadata is provided in the same form as supplied in E-MTAB-5522.
#' This contains information such as the cell type, plate of origin, spike-in addition order and oncogene induction. 
#'
#' Two sets of spike-ins were added to each cell in each dataset.
#' These are available as the \code{"SIRV"} and \code{"ERCC"} entries in the \code{\link{altExps}}.
#'
#' If \code{split.oncogene=TRUE} and \code{which="416b"},
#' the CBFB-MYH11-mcherry oncogene is moved to extra \code{"oncogene"} entry in the \code{\link{altExps}}.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
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
LunSpikeInData <- function(which=c("416b", "tropho"), split.oncogene=FALSE, location=TRUE, legacy=FALSE) {
    which <- match.arg(which)

    if (!legacy && !split.oncogene) {
        sce <- fetchDataset("lun-spikein-2017", "2023-12-18", path=which, realize.assays=TRUE) 

    } else {
        version <- "2.0.0"
        sce <- .create_sce(file.path("lun-spikein", version), suffix=which)

        spike.type <- rep("endogenous", nrow(sce))
        spike.type[grep("ERCC", rownames(sce))] <- "ERCC"
        spike.type[grep("SIRV", rownames(sce))] <- "SIRV"

        if (split.oncogene) {
            spike.type[rownames(sce)=="CBFB-MYH11-mcherry"] <- "oncogene"
        }

        sce <- splitAltExps(sce, spike.type, ref="endogenous")

        spike.exp <- altExp(sce, "ERCC")
        spikedata <- ERCCSpikeInConcentrations(volume = 100, dilution = 3e6)
        spikedata <- spikedata[rownames(spike.exp), ]
        rowData(spike.exp) <- cbind(rowData(spike.exp), spikedata)
        altExp(sce, "ERCC") <- spike.exp
    }

    .define_location_from_ensembl(sce, species="Mm", location=location)
}
