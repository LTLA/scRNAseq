# library(testthat); library(scRNAseq); source('test-pollen-glia.R')

test_that("Pollen glia getter works as expected", {
    sce <- PollenGliaData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- PollenGliaData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
