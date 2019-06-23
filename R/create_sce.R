#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.create_sce <- function(dataset, hub=ExperimentHub(), assay="counts", has.rowdata=TRUE, suffix=NULL) {
    host <- file.path("scRNAseq", dataset)
    if (is.null(suffix)) {
        suffix <- ""
    } else {
        suffix <- paste0("-", suffix)
    }

    all.assays <- list()
    for (a in assays) {
        all.assays[[a]] <- hub[hub$rdatapath==file.path(host, sprintf("%s%s.rds", a, suffix))][[1]]
    }

    coldata <- hub[hub$rdatapath==file.path(host, sprintf("coldata%s.rds", suffix))][[1]]
    if (has.rowdata) {
        rowdata <- hub[hub$rdatapath==file.path(host, sprintf("rowdata.rds", suffix))][[1]]
    } else {
        rowdata <- DataFrame(row.names=rownames(counts))
    }

    SingleCellExperiment(all.assays, rowData=rowdata, colData=coldata)
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
    SingleCellExperiment(all.assays, colData=coldata, metadata=metadata)
}
