write.csv(file="../extdata/metadata-leng-esc.csv", 
    data.frame(
        Title = sprintf("Leng ESC %s", c("normcounts", "colData")),
        Description = sprintf("%s for the Leng embryonic stem cell single-cell RNA-seq dataset", 
            c("Normalized count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "leng-esc", "2.0.0", c("normcounts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="hg19",
        SourceType="CSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64016",
        SourceVersion="GSE64016_H1andFUCCI_normalized_EC.csv.gz",
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE
)
