collected <- list()
for (analysis in c("Cellranger", "emptyDrops")) {
    collected[[analysis]] <- data.frame(
            Title = sprintf("Ernst spermatogenesis %s (%s cell calls)", c("counts", "colData", "rowData"), analysis),
            Description = sprintf("%s for Ernst mouse spermatogenesis single-cell RNA-seq dataset with the %s-derived cell calls", 
                c("Count matrix", "Per-cell metadata", "Per-gene metadata"), analysis),
            RDataPath = file.path("scRNAseq", "he-organ-atlas", "2.6.0", 
                paste0(c("counts", "coldata", "rowdata"), "-", tolower(analysis), ".rds")),
            BiocVersion="3.13",
            Genome="mm10",
            SourceType=c("MTX", "TXT", "TSV"),
            SourceUrl="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6946",
            SourceVersion=if (analysis=="Cellranger") {
                    c("E-MTAB-6946.processed.3.zip", 
                    "E-MTAB-6946.processed.1.zip", 
                    "E-MTAB-6946.processed.2.zip")
                } else {
                    c("E-MTAB-6946.processed.4.zip",
                    "E-MTAB-6946.processed.5.zip",
                    "E-MTAB-6946.processed.6.zip")
                },
            Species="Mus musculus",
            TaxonomyId="10090",
            Coordinate_1_based=NA,
            DataProvider="ArrayExpress",
            Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
            RDataClass=c("dgCMatrix", "DFrame", "DFrame"),
            DispatchClass="Rds",
            stringsAsFactors = FALSE
        )
}

write.csv(file="../../extdata/2.6.0/metadata-ernst-spermatogenesis.csv", rbind(collected[[1]], collected[[2]]), row.names=FALSE)
