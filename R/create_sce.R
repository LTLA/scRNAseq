#' @import methods
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors mcols
.get_hub_resource_by_rdatapath <- function(hub, rdata_dirname, rdata_basename) {
    stopifnot(is(hub, "ExperimentHub"),
              is.character(rdata_dirname),
              length(rdata_dirname) == 1L,
              !is.na(rdata_dirname),
              is.character(rdata_basename),
              length(rdata_basename) == 1L,
              !is.na(rdata_basename)
    )
    rdatapath <- file.path(rdata_dirname, rdata_basename)
    idx <- which(rdatapath == mcols(hub)$rdatapath)
    if (length(idx) == 0L)
        stop("rdatapath \"", rdatapath, "\" matches no ",
             class(hub), " resource")
    if (length(idx) > 1L)
        warning("rdatapath \"", rdatapath, "\" matches\n  more than one ",
                class(hub), " resource; picked the first one")
    updateObject(hub[[idx]])
}

#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce <- function(dataset, hub=.ExperimentHub(), assays="counts", has.rowdata=TRUE, has.coldata=TRUE, suffix=NULL) {
    host <- file.path("scRNAseq", dataset)
    if (is.null(suffix)) {
        suffix <- ""
    } else {
        suffix <- paste0("-", suffix)
    }

    all.assays <- list()
    for (a in assays) {
        rdata_basename <- sprintf("%s%s.rds", a, suffix)
        all.assays[[a]] <- .get_hub_resource_by_rdatapath(hub, host, rdata_basename)
    }

    args <- list()
    if (has.coldata) {
        rdata_basename <- sprintf("coldata%s.rds", suffix)
        args$colData <- .get_hub_resource_by_rdatapath(hub, host, rdata_basename)
    }
    if (has.rowdata) {
        rdata_basename <- sprintf("rowdata%s.rds", suffix)
        args$rowData <- .get_hub_resource_by_rdatapath(hub, host, rdata_basename)
    }

    do.call(SingleCellExperiment, c(list(assays=all.assays), args))
}

#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce_legacy <- function(dataset, assays, hub=.ExperimentHub()) {
    choices <- c("tophat_counts", "cufflinks_fpkm", "rsem_counts", "rsem_tpm")
    if (is.null(assays)) {
        assays <- choices
    } else {
        assays <- match.arg(assays, choices, several.ok=TRUE)
    }

    host <- file.path("scRNAseq", dataset)
    all.assays <- list()
    for (i in assays) {
        rdata_basename <- sprintf("%s.rds", i)
        all.assays[[i]] <- .get_hub_resource_by_rdatapath(hub, host, rdata_basename)
    }

    coldata <- .get_hub_resource_by_rdatapath(hub, host, "coldata.rds")
    metadata <- .get_hub_resource_by_rdatapath(hub, host, "metadata.rds")
    SingleCellExperiment(all.assays, colData=coldata, metadata=metadata)
}

.ExperimentHub <- function() {
    .move_cache("ExperimentHub", "EXPERIMENT_HUB_CACHE")
    ExperimentHub()
}
