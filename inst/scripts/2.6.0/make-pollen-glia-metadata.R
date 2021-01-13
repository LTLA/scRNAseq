write.csv(file="../../extdata/2.6.0/metadata-pollen-glia.csv",
    data.frame(
        Title = sprintf("Pollen radial glia %s", c("counts", "colData")),
        Description = sprintf("%s for the Pollen radial glia single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "pollen-glia", "2.6.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.13",
        Genome="hg19",
        SourceType=c("TSV", "XLSX"),
        SourceUrl="https://www.pollenlab.org/datasets",
        SourceVersion=c(
            "oRG paper - counts.txt",
            "Pollen%20et%20al%202015%20updated%20metadata.xlsx"),
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
