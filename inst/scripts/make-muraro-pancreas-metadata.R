write.csv(file="../extdata/metadata-muraro-pancreas.csv",
    data.frame(
        Title = sprintf("Muraro pancreas %s", c("counts", "rowData", "colData")),
        Description = sprintf("%s for the Muraro pancreas single-cell RNA-seq dataset", 
            c("Count matrix", "Per-gene metadata", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "muraro-pancreas", "2.0.0", 
            c("counts.rds", "rowdata.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="hg19",
        SourceType="TSV",
        SourceUrl=c(
            rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241", 2),
            "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/muraro"
        ),
        SourceVersion=c(
            rep("GSE85241_cellsystems_dataset_4donors_updated.csv.gz", 2),
            "cell_type_annotation_Cels2016.csv"
        ),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=TRUE,
        DataProvider="Mauro Javier Muraro",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="character",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
