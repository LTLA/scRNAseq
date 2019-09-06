write.csv(file="../extdata/metadata-nestorowa-hsc.csv", 
    data.frame(
        Title = sprintf("Nestorowa HSC %s", c("counts", "colData")),
        Description = sprintf("%s for the Nestorowa haematopoietic stem cell single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "nestorowa-hsc", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl=c(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682",
            "http://blood.stemcells.cam.ac.uk/data"
        ),
        SourceVersion=c(
            "GSE81682_HTSeq_counts.txt.gz",
            "all_cell_types.txt;coordinates_gene_counts_flow_cytometry.txt.gz"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider=c("GEO", "Bertie Göttgens"),
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
