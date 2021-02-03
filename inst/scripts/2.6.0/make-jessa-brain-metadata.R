write.csv(file="../../extdata/2.6.0/metadata-jessa-brain.csv",
    data.frame(
        Title = sprintf("Jessa mouse brain %s", c("counts", "colData", "rowData", "reducedDims")),
        Description = sprintf("%s for the Jessa mouse brain single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata", "Per-gene metadata", "Reduced dimensions")),
        RDataPath = file.path("scRNAseq", "jessa-brain", "2.6.0", 
            c("counts.rds", "coldata.rds", "rowdata.rds", "reddims.rds")),
        BiocVersion="3.13",
        Genome="mm10",
        SourceType=c("MTX", "TSV", "TSV", "TSV"),
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133531",
        SourceVersion=c(
            "GSE133531_RAW.tar",
            "GSE133531_RAW.tar;GSE133531_Forebrain_join.2D.tsv.gz;GSE133531_Pons_join.2D.tsv.gz"), # recycle
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame", "DFrame", "list"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
