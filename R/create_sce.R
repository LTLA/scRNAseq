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

#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce_legacy <- function(dataset, hub=ExperimentHub()) {
    host <- file.path("scRNAseq", dataset)

    assay.names <- c("tophat_counts", "cufflinks_fpkm", "rsem_counts", "rsem_tpm")
    all.assays <- list()
    for (i in assay.names) {
        path <- file.path(host, sprintf("%s.rds", i))
        all.assays[[i]] <- hub[hub$rdatapath==path][[1]]
    }

    coldata <- hub[hub$rdatapath==file.path(host, "coldata.rds")][[1]]
    metadata <- hub[hub$rdatapath==file.path(host, "metadata.rds")][[1]]
    SingleCellExperiment(all.assays, rowData=rowdata, colData=coldata, metadata=metadata)
}
