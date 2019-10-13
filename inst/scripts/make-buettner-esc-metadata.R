write.csv(file="../extdata/metadata-buettner-esc.csv", 
    data.frame(
        Title = sprintf("Buettner ESC %s", c("counts", "rowData", "colData")),
        Description = sprintf("%s for the Buettner embryonic stem cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-gene metadata", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "buettner-esc", "2.0.0", c("counts.rds", "rowdata.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="ZIP",
        SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805",
        SourceVersion="E-MTAB-2805.processed.1.zip",
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="ArrayExpress",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DFrame", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
