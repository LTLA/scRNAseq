GENERATOR <- function(short, long, target, src) {
    data.frame(
        Title = sprintf("Bach mammary %s", short),
        Description = sprintf("%s for the Bach mammary gland single-cell RNA-seq dataset", long), 
        RDataPath = file.path("scRNAseq", "bach-mammary", "2.0.0", target), 
        BiocVersion="3.10",
        Genome="mm10",
        SourceType=if (grepl("counts", target)) "MTX" else "TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106273",
        SourceVersion=src,
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass=if (grepl("counts", target)) "dgCMatrix" else "DFrame",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    )
}

accessions <- c(
    "GSM2834498", "GSM2834499",
    "GSM2834500", "GSM2834501",
    "GSM2834502", "GSM2834503",
    "GSM2834504", "GSM2834505")

samples <- c(
    "NP_1", "NP_2",
    "G_1", "G_2",
    "L_1", "L_2",
    "PI_1", "PI_2")

conditions <- c(NP="nulliparous",
    G="gestation",
    L="lactation",
    PI="post-involution")[sub("_.*", "", samples)]

collected <- list()
for (i in seq_along(accessions)) {
    n <- sub(".*_", "", samples[i])
    collected[[i]] <- rbind(
        GENERATOR(paste(samples[i], "counts"), paste0("Count matrix (", conditions[i], " replicate ", n, ")"), 
            paste0("counts-", samples[i], ".rds"), paste0(accessions[i], "_", samples[i], "_counts.mtx.gz")),
        GENERATOR(paste(samples[i], "colData"), paste0("Per-cell metadata (", conditions[i], " replicate ", n, ")"), 
            paste0("coldata-", samples[i], ".rds"), paste0(accessions[i], "_", samples[i], "_barcodes.tsv.gz"))
    )
} 
collected[[i+1]] <- GENERATOR("rowData", "Per-gene metadata", "rowdata.rds", 
    paste0(accessions[1], "_", samples[1], "_genes.tsv.gz"))

write.csv(file="../extdata/2.0.0/metadata-bach-mammary.csv", row.names=FALSE, do.call(rbind, collected))
