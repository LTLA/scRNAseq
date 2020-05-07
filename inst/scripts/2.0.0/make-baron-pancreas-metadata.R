FUN <- function(species, short, long, id) {
    data.frame(
        Title = sprintf("Baron %s pancreas %s", species, c("counts", "colData")),
        Description = sprintf("%s for the Baron %s pancreas single-cell RNA-seq dataset", 
            c("Count matrix", "Per-cell metadata"), species),
        RDataPath = file.path("scRNAseq", "baron-pancreas", "2.0.0", 
            sprintf(c("counts-%s.rds", "coldata-%s.rds"), species)),
        BiocVersion="3.10",
        Genome=short,
        SourceType="CSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133",
        SourceVersion="GSE84133_RAW.tar", 
        Species=long,
        TaxonomyId=id,
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=c("dgCMatrix", "DataFrame"),
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    )
}

write.csv(file="../extdata/2.0.0/metadata-baron-pancreas.csv", 
    rbind(FUN("human", "hg38", "Homo sapiens", "9606"), FUN("mouse", "mm10", "Mus musculus", "10090")),
    row.names=FALSE)
