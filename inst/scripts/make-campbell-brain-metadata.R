write.csv(file="../extdata/metadata-campbell-brain.csv", 
    data.frame(
        Title = sprintf("Campbell brain %s", c("counts", "colData")),
        Description = sprintf("%s for the Campbell brain single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "campbell-brain", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93374",
        SourceVersion=c(
            "GSE93374_Merged_all_020816_DGE.txt.gz",
            "GSE93374_cell_metadata.txt.gz"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
