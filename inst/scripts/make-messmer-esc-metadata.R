write.csv(file="../extdata/metadata-messmer-esc.csv", row.names=FALSE,
    data.frame(
        Title = sprintf("Messmer human ESC %s", c("counts", "rowData", "colData")),
        Description = sprintf("%s for the Messmer human embryonic stem cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-gene metadata", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "messmer-esc", "2.0.0", 
            c("counts.rds", "rowdata.rds", "coldata.rds")), 
        BiocVersion="3.10",
        Genome="hg38",
        SourceType="TSV",
        SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6819/",
        SourceVersion=c(
            "E-MTAB-6819.processed.1.zip",
            "E-MTAB-6819.processed.1.zip",
            "E-MTAB-6819.sdrf.txt"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="ArrayExpress",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DataFrame", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    )
)
