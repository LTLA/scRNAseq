write.csv(file="../extdata/metadata-paul-hsc.csv", 
    data.frame(
        Title = sprintf("Paul HSC %s", c("counts", "colData")),
        Description = sprintf("%s for the Paul haematopoietic stem cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "paul-hsc", "2.2.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.11",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72857",
        SourceVersion=c(
            "GSE72857_umitab.txt.gz",
            "GSE72857_experimental_design.txt.gz"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
