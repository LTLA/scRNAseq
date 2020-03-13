#' Obtain the ERCC concentration from Thermofisher
#'
#' @details
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ercc-spike-in-concentrations}.
#'
#' @param genes Character vector specifying the ERCC spike-in genes used.
#' @param volume Numeric scalar specifying the volume (in nanoliters) of ERCC 
#' spike-in mixture.
#' @param dilution Numeric scalar specifying the ratio of ERCC mixture to total 
#' well content.
#' @return A \code{data.frame} object with ERCC spike-in information.
#'
#' @author Alan O'Callaghan
#'
#' @examples
#' sce <- ERCCSpikeInConcentrations()
#' 
#' @export
#' @importFrom SingleCellExperiment splitAltExps
ERCCSpikeInConcentrations <- function(genes = NULL, volume = NULL, dilution = NULL) {
    version <- "2.0.0"
    hub <- ExperimentHub()
    host <- "scRNAseq/ercc-concentrations"
    file <- hub[hub$rdatapath==file.path(host, "2.2.0", "cms_095046.txt")][[1]]
    table <- read.delim(file, check.names = FALSE)
    if (!is.null(genes)) {
      ind <- match(genes, table[["ERCC ID"]])
      table <- table[ind, ]
    }
    if (!is.null(volume) & !is.null(dilution)) {    
      molarity <- table[["concentration in Mix 1 (attomoles/ul)"]] * (10^(-18))
      molecules <- molarity * (6.02214076 * (10^23))
      table$molarity_ul <- molarity
      table$molecules_ul <- molecules
      table$molecules <- table$molecules_ul / dilution
      table$molecules <- table$molecules * (volume / 1000)
    }
    table
}
