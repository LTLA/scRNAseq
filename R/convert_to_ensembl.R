#' @importFrom AnnotationDbi mapIds
#' @importFrom SummarizedExperiment rowData rowData<-
.convert_to_ensembl <- function(sce, symbols, species, ahub.id=NULL, keytype="SYMBOL", ensembl=FALSE, location=TRUE) {
    if (!ensembl) {
        return(sce)
    }

    edb <- .pull_down_ensdb(species, ahub.id=ahub.id)
    ensid <- mapIds(edb, keys=symbols, keytype=keytype, column="GENEID")
    keep <- !is.na(ensid) & !duplicated(ensid)
    if (!all(keep)) {
        sce <- sce[keep,]
    }

    old <- rownames(sce)
    rownames(sce) <- ensid[keep]
    rowData(sce)$originalName <- old

    .define_location_from_ensembl(sce, edb=edb, location=location)
}

#' @importFrom SummarizedExperiment rowRanges<- rowRanges rowData
#' @importFrom S4Vectors mcols<-
#' @importFrom BiocGenerics cbind
#' @importFrom methods as
#' @importClassesFrom GenomicRanges GRangesList
#' @importFrom GenomicFeatures genes
.define_location_from_ensembl <- function(sce, species, ahub.id=NULL, edb=NULL, location=TRUE) {
    if (!location) {
        return(sce)
    }

    if (is.null(edb)) {
        edb <- .pull_down_ensdb(species, ahub.id=ahub.id)
    }

    # Need to switch between a GRL and a GR, depending on whether
    # there are non-valid keys in 'sce'.
    ginfo <- genes(edb, columns=character(0))
    mcols(ginfo) <- NULL
    m <- match(rownames(sce), names(ginfo))
    present <- !is.na(m)

    if (!all(present)) {
        replacement <- rowRanges(sce)
        mcols(replacement) <- NULL
        ginfo <- as(ginfo, "GRangesList")
        replacement[present] <- ginfo[m[present]]
    } else {
        replacement <- ginfo[m]        
    }

    mcols(replacement) <- rowData(sce)
    rowRanges(sce) <- replacement
    sce
}

#' @importClassesFrom ensembldb EnsDb
#' @importFrom AnnotationHub AnnotationHub
.pull_down_ensdb <- function(species, ahub.id) {
    if (is.null(ahub.id)) {
        if (species=="Mm") {
            ahub.id <- "AH73905"
        } else if (species=="Hs") {
            ahub.id <- "AH73881"
        }
    }
    AnnotationHub()[[ahub.id]]
}
