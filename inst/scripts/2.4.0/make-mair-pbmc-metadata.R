write.csv(file="../../extdata/2.4.0/metadata-mair-pbmc.csv", 
    data.frame(
        Title = sprintf("Mair PBMC %s", c("RNA counts", "RNA rowData", "ADT counts", "ADT rowData", "colData")),
        Description = sprintf("%s for the Mair healthy PBMC targeted CITE-seq dataset", 
            c("RNA count matrix", "RNA gene metadata","ADT count matrix","ADT feature metadata", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "mair-pbmc", "2.4.0",
            c("counts-rna.rds", "rowdata-rna.rds", "counts-adt.rds", "rowdata-adt.rds","coldata.rds")),
        BiocVersion="3.12",
        Genome="hg38",
        SourceType="CSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135325",
        SourceVersion="May 19, 2020",
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Stephany Orjuela <sorjuelal@gmail.com>",
        RDataClass=c("matrix", "DFrame", "matrix", "DFrame", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
