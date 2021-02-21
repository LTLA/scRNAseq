write.csv(file="../../extdata/2.6.0/metadata-bunis-hspc.csv",
    data.frame(
        Title = sprintf("Bunis human HSPC %s", c("counts", "colData", "rowData")),
        Description = sprintf("%s for the Bunis human haematopoietic stem and progenitor single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata", "Per-gene metadata")),
        RDataPath = file.path("scRNAseq", "bunis-hspc", "2.6.0", 
            c("counts.rds", "coldata.rds", "rowdata.rds")),
        BiocVersion="3.13",
        Genome="hg38",
        SourceType=c("MTX", "TSV", "TSV"),
        SourceUrl=c(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086;https://ndownloader.figshare.com/files/25953740",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086"),
        SourceVersion=c(
            "GSE158490_matrix.mtx.gz",
            "GSE158490_HSPC.best.txt.gz;GSE158490_barcodes.tsv.gz",	
            "GSE158490_genes.tsv.gz"),
        Species="Homo sapiens",
        TaxonomyId="9606",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DFrame", "DFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
