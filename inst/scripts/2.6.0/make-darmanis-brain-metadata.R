write.csv(file="../../extdata/2.6.0/metadata-darmanis-brain.csv",
    data.frame(
        Title = sprintf("Darmanis brain %s", c("counts", "colData")),
        Description = sprintf("%s for the Darmanis human brain single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "darmanis-brain", "2.6.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.13",
        Genome="hg38",
        SourceType=c("TSV", "SOFT"),
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086",
        SourceVersion=c(
            "GSE67835_RAW.tar",
            "GSE67835-GPL15520_series_matrix.txt.gz;GSE67835-GPL18573_series_matrix.txt.gz"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
