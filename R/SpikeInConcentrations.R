#' Obtain the ERCC concentration from Thermofisher
#'
#' @details
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ercc-spike-in-concentrations}.
#'
#' @return A \code{data.frame} object with ERCC spike-in information.
#'
#' @author Alan O'Callaghan
#'
#' @examples
#' sce <- ERCCSpikeInConcentrations()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
ERCCSpikeInConcentrations <- function(volume_nl, dilution, genes, mix1 = TRUE) {
    version <- "2.0.0"
    hub <- ExperimentHub()
    host <- "scRNAseq/ercc-spike-in-concentrations"
    rdata <- hub[hub$rdatapath==file.path(host, "1.0.0", "concentrations.rds")][[1]]
    spike_table <- readRDS(rdata)

    molarity <- spike_table[["concentration in Mix 1 (attomoles/ul)"]] * (10^(-18))
    molecules <- molarity * (6.02214076 * (10^23))
    table$molarity_ul <- molarity
    table$molecules_ul <- molecules

    ind <- match(genes, spike_table[["ERCC.ID"]])
    spike_table <- spike_table[ind, ]

    count <- spike_table$molecules_ul / dilution
    count * (volume_nl / 1000)
}
