write.csv(file="../extdata/metadata-richard-tcell.csv", 
    data.frame(
        Title = sprintf("Richard mouse CD8+ T cell %s", c("counts", "colData")),
        Description = sprintf("%s for the Richard mouse CD8+ T cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "richard-tcell", "2.0.0", c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6051/",
        SourceVersion=c("E-MTAB-6051.processed.1.zip", "E-MTAB-6051.sdrf.txt"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="ArrayExpress",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
