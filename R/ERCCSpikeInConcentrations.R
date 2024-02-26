#' Obtain ERCC concentrations
#'
#' Obtain ERCC spike-in concentrations from the Thermo Fisher Scientific website.
#'
#' @param volume Numeric scalar specifying the added volume (in nanoliters) of ERCC spike-in mixture.
#' Only used if \code{dilution} is specified.
#' @param dilution Numeric scalar specifying the dilution factor used for the added volume of the spike-in mixture.
#' Only used if \code{volume} is specified.
#' @param mix String specifying whether to compute the number of molecules for mix 1 or 2.
#' Only used if both \code{dilution} and \code{volume} are specified.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' If \code{volume} and \code{dilution} are specified,
#' an additional column is added to the output specifying the number of molecules of spike-in transcipt for the specified mix.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/ercc-concentrations}.
#'
#' @return A \linkS4class{DataFrame} object with one row per ERCC spike-in transcript.
#' This contains information such as the spike-in concentration in each mix.
#'
#' @author Alan O'Callaghan
#'
#' @examples
#' df <- ERCCSpikeInConcentrations()
#' 
#' @export
#' @importFrom utils read.delim
#' @importFrom SingleCellExperiment splitAltExps
ERCCSpikeInConcentrations <- function(volume = NULL, dilution = NULL, mix=c("1", "2"), legacy=FALSE) {
    if (!legacy) {
        table <- fetchDataset("ercc", "2023-12-20")
    } else {
        version <- "2.0.0"
        hub <- .ExperimentHub()
        host <- "scRNAseq/ercc-concentrations"
        file <- hub[hub$rdatapath==file.path(host, "2.2.0", "cms_095046.txt")][[1]]

        table <- read.delim(file, check.names = FALSE)
        rownames(table) <- table[,"ERCC ID"]
        table <- table[,-(1:2)]
    }

    if (!is.null(volume) && !is.null(dilution)) {
        field <- sprintf("concentration in Mix %s (attomoles/ul)", match.arg(mix))
        molarity <- table[[field]] * 10^-18 # get concentration in mol/uL.
        molecules_ul <- molarity * (6.02214076 * (10^23))
        table$molecules <- molecules_ul / dilution * (volume / 1000) # get volume in uL.
    }

    DataFrame(table, check.names=FALSE)
}
