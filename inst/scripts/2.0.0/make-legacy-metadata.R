base <- c("Tophat count", "Cufflinks FPKM", "RSEM count", "RSEM TPM")
file.titles <- paste0(base, "s")
file.descs <- paste(base, "matrix")
file.paths <- paste0(c("tophat_counts", "cufflinks_fpkm", "rsem_counts", "rsem_tpm"), ".rds")

CREATOR <- function(desc, path, genome, src, species, taxid, append=TRUE) {
    write.table(file="../extdata/2.0.0/metadata-legacy.csv",
        data.frame(
            Title = sprintf("%s %s", desc, c(file.titles, "colData", "metadata")),
            Description = sprintf("%s for the %s single-cell RNA-seq dataset", 
                c(file.descs, "Per-cell metadata", "Experiment metadata"), desc),
            RDataPath = file.path("scRNAseq", path, c(file.paths, "coldata.rds", "metadata.rds")),
            BiocVersion="3.10",
            Genome=genome,
            SourceType="RData",
            SourceUrl="http://bioconductor.org/packages/3.9/data/experiment/src/contrib/scRNAseq_1.10.0.tar.gz",
            SourceVersion=src,
            Species=species,
            TaxonomyId=taxid,
            Coordinate_1_based=NA,
            DataProvider="Michael Cole",
            Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
            RDataClass=c(rep("matrix", length(file.titles)), "DataFrame", "list"),
            DispatchClass="Rds",
            stringsAsFactors = FALSE
        ),
        row.names=FALSE, sep=",", append=append, col.names=!append)
}

CREATOR("Allen brain", "legacy-allen/1.10.0", "mm10", "data/allen.rda", "Mus musculus", "10090", append=FALSE)
CREATOR("T helper 2", "legacy-th2/1.10.0", "mm10", "data/th2.rda", "Mus musculus", "10090")
CREATOR("Fluidigm", "legacy-fluidigm/1.10.0", "hg38", "data/pollen.rda", "Homo sapiens", "9606")
