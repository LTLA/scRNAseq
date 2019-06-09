#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.create_sce <- function(dataset, hub=ExperimentHub(), has.rowdata=TRUE) {
    host <- file.path("scRNAseq", dataset)
    counts <- hub[hub$rdatapath==file.path(host, "counts.rds")][[1]]
    coldata <- hub[hub$rdatapath==file.path(host, "coldata.rds")][[1]]
    if (has.rowdata) {
        rowdata <- hub[hub$rdatapath==file.path(host, "rowdata.rds")][[1]]
    } else {
        rowdata <- DataFrame(row.names=rownames(counts))
    }
    SingleCellExperiment(list(counts=counts), rowData=rowdata, colData=coldata)
}
