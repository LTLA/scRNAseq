#' Obtain the ERCC concentration from Thermofisher
#'
#' @details
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ercc-spike-in-concentrations}.
#'
#' @param genes A character vector specifying the ERCC spike-in genes used.
#' @param volume_nl Volume (in nanoliters) of ERCC spike-in mixture.
#' @param dilution The ratio of ERCC mixture to total well content.
#' @return A \code{data.frame} object with ERCC spike-in information.
#'
#' @author Alan O'Callaghan
#'
#' @examples
#' sce <- ERCCSpikeInConcentrations()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
ERCCSpikeInConcentrations <- function(genes = NULL, volume_nl = NULL, dilution = NULL) {
    version <- "2.0.0"
    hub <- ExperimentHub()
    host <- "scRNAseq/ercc-spike-in-concentrations"
    rdata <- hub[hub$rdatapath==file.path(host, "1.0.0", "concentrations.rds")][[1]]
    out <- readRDS(rdata)

    if (!is.null(genes)) {
      ind <- match(genes, spike_table[["ERCC.ID"]])
      out <- out[ind, ]
    }
    if (!is.null(volume_nl) & !is.null(dilution)) {    
      molarity <- spike_table[["concentration in Mix 1 (attomoles/ul)"]] * (10^(-18))
      molecules <- molarity * (6.02214076 * (10^23))
      out$molarity_ul <- molarity
      out$molecules_ul <- molecules
      out$molecules <- out$molecules_ul / dilution
      out$molecules <- out$molecules * (volume_nl / 1000)
    }
    out
}
