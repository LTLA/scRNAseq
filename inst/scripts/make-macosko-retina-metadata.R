write.csv(file="../extdata/metadata-macosko-retina.csv", 
    data.frame(
        Title = sprintf("Macosko retina %s", c("counts", "colData")),
        Description = sprintf("%s for the Macosko retina single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "macosko-retina", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl=c(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472",
            "http://mccarrolllab.com/wp-content/uploads/2015/05"
        ),
        SourceVersion=c(
            "GSE63472_P14Retina_merged_digital_expression.txt.gz",
            "retina_clusteridentities.txt"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider=c("GEO", "Steve McCarroll"),
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
