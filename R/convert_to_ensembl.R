.convert_to_ensembl <- function(sce, symbols, species, keytype="SYMBOL") {
    if (species=="Mm") {
        tag <- "AH73905"
    } else if (species=="Hs") {
        tag <- "AH73881"
    }
    edb <- AnnotationHub::AnnotationHub()[[tag]]
    ensid <- AnnotationDbi::mapIds(edb, keys=symbols, keytype=keytype, column="GENEID")

    keep <- !is.na(ensid) & !duplicated(ensid)
    if (!all(keep)) {
        sce <- sce[keep,]
    }
    rownames(sce) <- ensid[keep]
    sce
}
