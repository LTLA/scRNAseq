# library(testthat); library(scRNAseq); source('test-he-organ-atlas.R')

test_that("He organ atlas getter works as expected", {
    sce <- HeOrganAtlasData()
    expect_s4_class(sce, "SingleCellExperiment")

    suppressWarnings(sce2 <- HeOrganAtlasData(ensembl=TRUE))
    expect_true(nrow(sce2) < nrow(sce))

    solo.sce <- HeOrganAtlasData(tissue="Blood")
    expect_identical(reducedDimNames(solo.sce), "tSNE")
})
