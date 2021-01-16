write.csv(file="../../extdata/2.6.0/metadata-bhaduri-organoid.csv",
    data.frame(
        Title = sprintf("Bhaduri organoid %s", c("counts", "colData")),
        Description = sprintf("%s for the Bhaduri cortical organoid single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "bhaduri-organoid", "2.6.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.13",
        Genome="hg38",
        SourceType=c("TSV", "TXT"),
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132672",
        SourceVersion=c(
            "GSE132672_allorganoids_withnew_matrix_txt.gz",
            "GSE132672_series_matrix.txt.gz"),
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
