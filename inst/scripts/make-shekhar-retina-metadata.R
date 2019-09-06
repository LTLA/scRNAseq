write.csv(file="../extdata/metadata-shekhar-retina.csv", 
    data.frame(
        Title = sprintf("Shekhar retina %s", c("counts", "colData")),
        Description = sprintf("%s for the Shekhar retina single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata")),
        RDataPath = file.path("scRNAseq", "shekhar-retina", "2.0.0", 
            c("counts.rds", "coldata.rds")),
        BiocVersion="3.10",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl=c(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81904",
            "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/shekhar"
        ),
        SourceVersion=c(
            "GSE81904_BipolarUMICounts_Cell2016.txt.gz",
            "clust_retinal_bipolar.txt"),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider=c("GEO", "Martin Hemberg"),
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("matrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
