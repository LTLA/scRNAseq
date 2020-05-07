write.csv(file="../extdata/2.0.0/metadata-lawlor-pancreas.csv",
    data.frame(
        Title = sprintf("Lawlor pancreas %s", c("counts", "colData")),
        Description = sprintf("%s for the Lawlor pancreas single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "lawlor-pancreas", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="hg19",
        SourceType="TXT",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469",
        SourceVersion=c(
            "GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz",
            "GSE86469_series_matrix.txt.gz"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
