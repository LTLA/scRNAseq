write.csv(file="../../extdata/2.8.0/metadata-ledergor-myeloma.csv",
    data.frame(
        Title = sprintf("Ledergor multiple myeloma %s", c("counts", "colData")),
        Description = sprintf("%s for the Ledergor multiple myeloma single-cell RNA-seq dataset",
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "ledergor-myeloma", "2.8.0",
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.14",
        Genome="hg38",
        SourceType=c("MTX", "TSV"),
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117156",
        SourceVersion=c(
            "GSE117156_RAW.tar",
            "GSE117156_metadata.txt.gz"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
