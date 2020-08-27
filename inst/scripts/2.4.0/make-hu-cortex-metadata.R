samples <- c(
    "cell-3T3",
    "nuclei-3T3",
    "nuclei-ctx-1",
    "nuclei-ctx-2",
    "nuclei-ctx-3", 
    "nuclei-ctx-4",
    "nuclei-ctx-5",
    "nuclei-ctx-6",
    "nuclei-ctx-7",
    "nuclei-ctx-8",
    "nuclei-ctx-9",
    "nuclei-ctx-10",
    "nuclei-ctx-11",
    "nuclei-ctx-12",
    "nuclei-ctx-13",
    "nuclei-ctx-saline1", 
    "nuclei-ctx-PTZ1", 
    "nuclei-ctx-saline2",
    "nuclei-ctx-PTZ2" 
)

details <- c(
    "3T3 cells",
    "3T3 nuclei",
    sprintf("mouse cortex cells (animal %i)", seq_len(13)),
    sprintf("mouse cortex cells (%s-treated, animal %i)",
        rep(c("saline", "PTZ"), rep=2),
        rep(1:2, each=2)
    )
)

files <- paste0("GSM28453", 58 + seq_along(samples))

write.csv(file="../../extdata/2.4.0/metadata-hu-cortex.csv", 
    data.frame(
        Title = sprintf("Hu %s counts", sub("cortex cells", "cortex", details)),
        Description = sprintf("Count matrix for %s in the Hu single-nucleus RNA-seq dataset", details),
        RDataPath = file.path("scRNAseq", "hu-cortex", "2.4.0", sprintf("counts-%s.rds", samples)),
        BiocVersion="3.12",
        Genome="mm10",
        SourceType="TSV",
        SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106678",
        SourceVersion=sprintf("%s_%s_gene_exonic.intronic_tagged.dge.txt.gz", files, samples),
        Species="Mus musculus",
        TaxonomyId="10090",
        Coordinate_1_based=NA,
        DataProvider="GEO",
        Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
        RDataClass="dgCMatrix",
        DispatchClass="Rds",
        stringsAsFactors = FALSE
    ),
    row.names=FALSE)
