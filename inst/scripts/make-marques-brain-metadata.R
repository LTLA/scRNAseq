write.csv(file="../extdata/metadata-marques-brain.csv", 
    data.frame(
        Title = sprintf("Marques brain %s", c("counts", "colData")),
        Description = sprintf("%s for the Marques brain single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "marques-brain", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75330",
        SourceVersion=c(
            "GSE75330_Marques_et_al_mol_counts2.tab.gz",
            "GSE75330_series_matrix.txt.gz"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
