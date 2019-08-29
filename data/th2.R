.Deprecated(msg="'data(th2)' is deprecated.\nUse ReprocessedTh2Data() instead.")
th2 <- (function(){
    bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
    bpath <- BiocFileCache::bfcrpath(bfc, "https://github.com/LTLA/scRNAseq/raw/80b3d4dde812b72406499922c5cd157e7b8d3b45/data/th2.rda")
    env <- new.env()
    load(bpath, envir=env)
    env$th2
})()
