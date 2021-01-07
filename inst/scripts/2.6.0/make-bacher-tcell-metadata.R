write.csv(file="../../extdata/2.6.0/metadata-bacher-tcell.csv",
    data.frame(
        Title = sprintf("Bacher COVID T cell %s", c("counts", "colData")),
        Description = sprintf("%s for the Bacher COVID T cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "bacher-tcell", "2.6.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.13",
        Genome="hg38",
        SourceType="TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086",
        SourceVersion=c(
            "GSE162086_RAW.tar",
            "GSE162086_seurat_metadata.tsv.gz"),
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
