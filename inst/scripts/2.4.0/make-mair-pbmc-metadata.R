write.csv(file="../extdata/metadata-mair-pbmc.csv", 
    data.frame(
        Title = sprintf("Mair PBMCs %s", c("counts", "colData")),
        Description = sprintf("%s for the Mair healthy PBMCs targetted CITEseq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("CITEseq", "mair-pbmc", c("counts.rds", "coldata.rds")),
        BiocVersion="3.11",
        Genome="hg38",
        SourceType="ZIP",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135325",
        SourceVersion="May 19, 2020",
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="ArrayExpress",
        Maintainer="Stephany Orjuela <sorjuelal@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
