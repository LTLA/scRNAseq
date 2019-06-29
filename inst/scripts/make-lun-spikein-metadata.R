FUN <- function(type, short) {
    data.frame(
        Title = sprintf("Lun %s plus spike-in %s", type, c("counts", "rowData", "colData")),
        Description = sprintf("%s for the Lun %s (plus spike-ins) single-cell RNA-seq dataset", 
            c("Count matrix", "Per-gene metadata", "Per-cell metadata"), type),
        RDataPath = file.path("scRNAseq", "lun-spikein", "2.0.0", 
            sprintf(c("counts-%s.rds", "rowdata-%s.rds", "coldata-%s.rds"), short)),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5522/",
        SourceVersion=c(
            "E-MTAB-5522.processed.1.zip",
            "E-MTAB-5522.processed.1.zip",
            "E-MTAB-5522.sdrf.txt"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=TRUE,
        DataProvider="Aaron Lun",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="character",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    )
}

write.csv(file="../extdata/metadata-lun-spikein.csv", 
    rbind(FUN("416B", "416b"), FUN("trophoblast", "tropho")),
    row.names=FALSE)
