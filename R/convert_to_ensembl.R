.convert_to_ensembl <- function(sce, symbols, species, keytype="SYMBOL", ensembl=FALSE, location=TRUE) {
    if (!ensembl) {
        return(sce)
    }

    edb <- .pull_down_ensdb(species)

    ensid <- AnnotationDbi::mapIds(edb, keys=symbols, keytype=keytype, column="GENEID")
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
.define_location_from_ensembl <- function(sce, species, edb=NULL, location=TRUE) {
    if (!location) {
        return(sce)
    }

    if (is.null(edb)) {
        edb <- .pull_down_ensdb(species)
    }

    # Need to switch between a GRL and a GR, depending on whether
    # there are non-valid keys in 'sce'.
    ginfo <- GenomicFeatures::genes(edb, columns=character(0))
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

.pull_down_ensdb <- function(species) {
    if (species=="Mm") {
        tag <- "AH73905"
    } else if (species=="Hs") {
        tag <- "AH73881"
    }
    AnnotationHub::AnnotationHub()[[tag]]
}

