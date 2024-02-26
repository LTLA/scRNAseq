#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce <- function(dataset, hub=ExperimentHub(), assays="counts", has.rowdata=TRUE, has.coldata=TRUE, suffix=NULL) {
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

    args <- list()
    if (has.coldata) {
        args$colData <- hub[hub$rdatapath==file.path(host, sprintf("coldata%s.rds", suffix))][[1]]
    }
    if (has.rowdata) {
        args$rowData <- hub[hub$rdatapath==file.path(host, sprintf("rowdata%s.rds", suffix))][[1]]
    }

    do.call(SingleCellExperiment, c(list(assays=all.assays), args))
}

#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce_legacy <- function(dataset, assays, hub=ExperimentHub()) {
    choices <- c("tophat_counts", "cufflinks_fpkm", "rsem_counts", "rsem_tpm")
    if (is.null(assays)) {
        assays <- choices
    } else {
        assays <- match.arg(assays, choices, several.ok=TRUE)
    }

    host <- file.path("scRNAseq", dataset)
    all.assays <- list()
    for (i in assays) {
        path <- file.path(host, sprintf("%s.rds", i))
        all.assays[[i]] <- hub[hub$rdatapath==path][[1]]
    }

    coldata <- hub[hub$rdatapath==file.path(host, "coldata.rds")][[1]]
    metadata <- hub[hub$rdatapath==file.path(host, "metadata.rds")][[1]]
    SingleCellExperiment(all.assays, colData=coldata, metadata=metadata)
}
