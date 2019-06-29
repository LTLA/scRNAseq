spawn_df <- function(type, suffix, src, file="TSV") {
    data.frame(
        Title = sprintf("La Manno %s %s", type, c("counts", "colData")),
        Description = sprintf("%s for the La Manno %s single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata"), type),
        RDataPath = file.path("scRNAseq", "lamanno-brain", "2.0.0", 
            sprintf(c("counts-%s.rds", "coldata-%s.rds"), suffix)),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType=file,
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76381",
        SourceVersion=src,
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=TRUE,
        DataProvider="Sten Linnarsson",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="character",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    )
}

write.csv(file="../extdata/metadata-lamanno-brain.csv", 
    rbind(
        spawn_df("human ES", "human-es", "GSE76381_ESMoleculeCounts.cef.txt.gz"),
        spawn_df("human embryo midbrain", "human-embryo", "GSE76381_EmbryoMoleculeCounts.cef.txt.gz"),
        spawn_df("mouse adult dopaminergic neuron", "mouse-adult", "GSE76381_MouseAdultDAMoleculeCounts.cef.txt.gz", file="CSV"),
        spawn_df("mouse embryo midbrain", "mouse-embryo", "GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz"),
        spawn_df("human iPS", "human-ips", "GSE76381_iPSMoleculeCounts.cef.txt.gz")
    ),
    row.names=FALSE)
