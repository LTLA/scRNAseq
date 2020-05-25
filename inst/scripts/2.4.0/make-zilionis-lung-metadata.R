FUN <- function(species, short, long, id) {
  data.frame(
    Title = sprintf("Zilionis %s lung %s", species, c("counts", "colData")),
    Description = sprintf("%s for the Zilionis %s lung single-cell RNA-seq dataset", 
                          c("Count matrix", "Per-cell metadata"), species),
    RDataPath = file.path("scRNAseq", "zilionis-lung", "2.4.0", 
                          sprintf(c("counts-%s.rds", "coldata-%s.rds"), species)),
    BiocVersion="3.12",
    Genome=short,
    SourceType="CSV",
    SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465",
    SourceVersion="GSE127465_RAW.tar", 
    Species=long,
    TaxonomyId=id,
    Coordinate_1_based=NA,
    DataProvider="GEO",
    Maintainer="Jens Preussner <jens.preussner@mpi-bn.mpg.de>",
    RDataClass=c("dgCMatrix", "DFrame"),
    DispatchClass="Rds",
    stringsAsFactors = FALSE
  )
}

write.csv(file="../../extdata/2.4.0/metadata-zilionis-lung.csv", 
    rbind(FUN("human", "hg38", "Homo sapiens", "9606"), FUN("mouse", "mm10", "Mus musculus", "10090")),
    row.names=FALSE)
