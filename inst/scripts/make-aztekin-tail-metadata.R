write.csv(file="../extdata/metadata-aztekin-tail.csv", 
    data.frame(
        Title = sprintf("Aztekin Xenopus tail %s", c("counts", "colData")),
        Description = sprintf("%s for the Aztekin Xenopus tail single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "aztekin-tail", "2.0.0", c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="v9.1",
        SourceType="ZIP",
        SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7716",
        SourceVersion="E-MTAB-7716.processed.1.zip",
        Species="Xenopus laevis",
        TaxonomyId="8355",
        Coordinate_1_based=NA,
        DataProvider="ArrayExpress",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
