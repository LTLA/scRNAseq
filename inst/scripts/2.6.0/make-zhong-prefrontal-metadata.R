write.csv(file="../../extdata/2.6.0/metadata-zhong-prefrontal.csv",
    data.frame(
        Title = sprintf("Zhong prefrontal cortex %s", c("counts", "colData")),
        Description = sprintf("%s for the Zhong prefrontal cortex single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "zhong-prefrontal", "2.6.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.13",
        Genome="hg19",
        SourceType=c("TSV", "XLSX"),
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276",
        SourceVersion=c(
            "GSE104276_all_pfc_2394_UMI_count_NOERCC",
            "GSE104276_readme_sample_barcode.xlsx	"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
