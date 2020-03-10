write.csv(file="../extdata/2.0.0/metadata-kolodziejczyk-esc.csv", 
    data.frame(
        Title = "Kolodziejczyk ESC counts",
        Description = "Count matrix for the Kolodziejczyk embryonic stem cell single-cell RNA-seq dataset", 
        RDataPath = file.path("scRNAseq", "kolodziejczyk-esc", "2.0.0", "counts.rds"),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TXT",
        SourceUrl="https://espresso.teichlab.sanger.ac.uk/static",
        SourceVersion="counttable_es.csv",
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="Sarah Teichmann",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="matrix", 
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
