write.csv(file="../../extdata/2.6.0/metadata-giladi-hsc.csv",
    data.frame(
        Title = sprintf("Giladi HSC %s", c("RNA counts", "RNA colData", "CRISPR counts", "CRISPR colData", "CRISPR rowData")),
        Description = sprintf("%s for the Giladi haematopoietic stem cell single-cell RNA-seq and CRISPR-seq dataset", 
            c("RNA count matrix", "Per-cell metadata (RNA)", "CRISPR count matrix", "Per-cell metadata (CRISPR)", "Per-guide metadata")),
        RDataPath = file.path("scRNAseq", "giladi-hsc", "2.6.0", 
            c("counts-rna.rds", "coldata-rna.rds", "counts-crispr.rds", "coldata-crispr.rds", "rowdata-crispr.rds")),
        BiocVersion="3.13",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl=rep(c(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113494",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113494"),
            c(2,3)),
        SourceVersion=c(
            "GSE92575_RAW.tar",
            "GSE92575_metadata.txt.gz",
            rep("GSE113494_crispseq_count.txt.gz", 3)),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame", "dgCMatrix", "DFrame", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
