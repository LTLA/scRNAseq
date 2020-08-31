write.csv(file="../../extdata/2.4.0/metadata-kotliarov-pbmc.csv", 
    data.frame(
        Title = sprintf("Kotliarov PBMCs %s", c("RNA counts", "ADT counts", "colData")),
        Description = sprintf("%s for the Kotliarov influenza-vaccinated healthy PBMCs CITE-seq dataset", 
            c("RNA count matrix", "ADT count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "kotliarov-pbmc", "2.4.0",
            c("counts-rna.rds", "counts-adt.rds", "coldata.rds")),
        BiocVersion="3.12",
        Genome="hg19",
        SourceType="RDS",
        SourceUrl="https://nih.figshare.com/ndownloader/files/20706642",
        SourceVersion="H1_day0_demultilexed_singlets.RDS",
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="NIH",
        Maintainer="Stephany Orjuela <sorjuelal@gmail.com>",
        RDataClass=c("dgCMatrix", "dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
