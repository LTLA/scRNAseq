CREATOR <- function(desc, path, genome, src, species, taxid, append=TRUE) {
    write.table(file="../extdata/metadata-legacy.csv",
        data.frame(
            Title = sprintf("%s %s", desc, c("counts", "colData", "metadata")),
            Description = sprintf("%s for the %s single-cell RNA-seq dataset", 
                c("Count matrix", "Per-cell metadata", "Experiment metadata"), desc),
            RDataPath = file.path("scRNAseq", path, c("counts.rds", "coldata.rds", "metadata.rds")),
            BiocVersion="3.10",
            Genome=genome,
            SourceType="TXT",
            SourceUrl="http://bioconductor.org/packages/3.9/data/experiment/src/contrib/scRNAseq_1.10.0.tar.gz",
            SourceVersion=src,
            Species=species,
            TaxonomyId=taxid,
            Coordinate_1_based=TRUE,
            DataProvider="Michael Cole",
            Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
            RDataClass="character",
            DispatchClass="Rds",
            stringsAsFactors = FALSE
        ),
        row.names=FALSE, sep=",", append=append, col.names=!append)
}

CREATOR("Allen brain", "legacy-allen/1.10.0", "mm10", "data/allen.rda", "Mus musculus", "10090", append=FALSE)
CREATOR("T helper 2", "legacy-th2/1.10.0", "mm10", "data/th2.rda", "Mus musculus", "10090")
CREATOR("Fluidigm", "legacy-fluidigm/1.10.0", "hg38", "data/pollen.rda", "Homo sapiens", "9606")
