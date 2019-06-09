write.csv(file="../extdata/metadata-tasic-brain-metadata.csv", 
    data.frame(
        Title = sprintf("Tasic brain %s", c("counts", "colData")),
        Description = sprintf("%s for the Tasic brain single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "tasic-brain", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TXT",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",
        SourceVersion=c(
            paste(
                "GSE71585_RefSeq_counts.csv.gz",
                "GSE71585_ERCC_and_tdTomato_counts.csv.gz", sep=";"),
            "GSE71585_Clustering_Results.csv.gz"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=TRUE,
        DataProvider="Lucas T Gray",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="character",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
