#' @importClassesFrom Matrix dgCMatrix
#' @importFrom SummarizedExperiment assay<- assay
realize_assays <- function(x) {
    for (y in assayNames(x)) {
        assay(x, y, withDimnames=FALSE) <- as(assay(x, y, withDimnames=FALSE), "dgCMatrix")
    }

    if (is(x, "SingleCellExperiment")) {
        for (z in altExpNames(x)) {
            altExp(x, z, withDimnames=FALSE) <- realize_assays(altExp(x, z, withDimnames=FALSE))
        }
    }

    x
}
