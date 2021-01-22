collected <- list()
for (tissue in c("Bladder", "Blood", "Common.bile.duct", "Esophagus", "Heart", 
        "Liver", "Lymph.node", "Marrow", "Muscle", "Rectum", "Skin",
        "Small.intestine", "Spleen", "Stomach", "Trachea")) {    
    nice.tissue <- tolower(sub("\\.", " ", tissue))
    collected[[tissue]] <- data.frame(
            Title = sprintf("He human organ atlas %s %s", nice.tissue, c("counts", "colData")),
            Description = sprintf("%s for the %s subset of the He human organ atlas single-cell RNA-seq dataset", 
                c("Count matrix", "Per-cell metadata"), nice.tissue),
            RDataPath = file.path("scRNAseq", "he-organ-atlas", "2.6.0", 
                paste0(c("counts", "coldata"), "-", tissue, ".rds")),
            BiocVersion="3.13",
            Genome="hg38",
            SourceType=c("CSV", "TSV"),
            SourceUrl=c(
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159929",
                "https://github.com/bei-lab/scRNA-AHCA/tree/master/Cell_barcode_and_corresponding_cell_types_of_AHCA"),
            SourceVersion=c(
                "GSE159929_RAW.tar",
                "Annotation_AHCA_alltissues_meta.data_84363_cell.txt;B_and_plasma.meta.data.txt;CD4_meta.data.txt;CD8_meta.data.txt;Endothelial_cell.meta.data.txt;Epithelial_cells.meta.data.txt;FibSmo.meta.data.txt;Myeloid.meta.data.txt"),
            Species="Homo sapiens",
            TaxonomyId="9606",
            Coordinate_1_based=NA,
            DataProvider="GEO",
            Maintainer="Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>",
            RDataClass=c("dgCMatrix", "DFrame"),
            DispatchClass="Rds",
            stringsAsFactors = FALSE
        )
}

write.csv(file="../../extdata/2.6.0/metadata-he-organ-atlas.csv", do.call(rbind, collected), row.names=FALSE)

