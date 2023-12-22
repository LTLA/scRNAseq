# This checks that saveDataset works as expected.
# library(testthat); library(scRNAseq); source("test-saveDataset.R")

sce <- SingleCellExperiment(list(counts=matrix(rpois(1000, lambda=1), 100, 10)))
sce$foo <- sample(letters, 10)
rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
colnames(sce) <- head(LETTERS, 10)

test_that("saveDataset works as expected", {
    meta <- list(
        title="My dataset",
        description="This is my dataset",
        taxonomy_id="10090",
        genome="GRCh38",
        sources=list(list(provider="GEO", id="GSE12345")),
        maintainer_name="Shizuka Mogami",
        maintainer_email="mogami.shizuka@765pro.com"
    )

    tmp <- tempfile()
    saveDataset(sce, tmp, meta)
    roundtrip <- alabaster.base::readObject(tmp)
    expect_identical(colData(roundtrip), colData(sce))
    expect_identical(as.matrix(counts(roundtrip)), assay(sce))

    meta <- jsonlite::fromJSON(file.path(tmp, "_bioconductor.json"), simplifyVector=FALSE)
    expect_identical(meta$bioconductor_version, as.character(BiocManager::version()))

    # Validation fails as expected.
    tmp <- tempfile()
    meta$title <- 1234
    expect_error(saveDataset(sce, tmp, meta), "title")
})
