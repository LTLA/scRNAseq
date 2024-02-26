#' Obtain ERCC molecule counts from concentrations
#'
#' Compute the number of molecules for each transcript in the ERCC spike-in mixture,
#' based on their published concentration as well as the volume of the diluted mixture added to each cell.
#'
#' @param volume Numeric scalar specifying the added volume (in microliters) of ERCC spike-in mixture.
#' @param dilution Numeric scalar specifying the dilution factor used for the added volume of the spike-in mixture.
#' @param mix String specifying whether to compute the number of molecules for mix 1 or 2.
#' @param ... Further arguments to pass to \code{\link{fetchDataset}}.
#'
#' @return A \linkS4class{DataFrame} object with one row per ERCC spike-in transcript.
#' This contains the estimated concentration and molecule count for each transcript.
#'
#' @author Aaron Lun, based on code from Alan O'Callaghan
#'
#' @examples
#' countErccMolecules(volume = 9, dilution = 300000)
#' 
#' @export
#' @importFrom S4Vectors DataFrame
countErccMolecules <- function(volume, dilution, mix=c("1", "2"), ...) {
    tab <- fetchDataset("ercc", "2023-12-20", ...)

    field <- sprintf("concentration in Mix %s (attomoles/ul)", match.arg(mix))
    molarity <- tab[[field]] * 10^-18 # get concentration in mol/uL.
    molecules_ul <- molarity * (6.02214076 * (10^23))

    # Previously, 'volume' was treated as nanolitres, but I'm pretty sure that
    # was not correct; who pipettes nL anyway? I'm pretty sure uL was intended.
    molecules <- molecules_ul / dilution * volume

    DataFrame(
        row.names=rownames(tab),
        subgroup=tab$subgroup,
        concentration=tab[[field]],
        molecules=molecules
    )
}
