.Deprecated(msg="'data(fluidigm)' is deprecated.\nUse ReprocessedFluidigmData() instead.")
bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
bpath <- BiocFileCache::bfcrpath(bfc, "https://github.com/LTLA/scRNAseq/raw/80b3d4dde812b72406499922c5cd157e7b8d3b45/data/fluidigm.rda")
env <- new.env()
load(bpath, envir=env)
fluidigm <- env$fluidigm
