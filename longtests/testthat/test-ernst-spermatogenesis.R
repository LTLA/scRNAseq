# library(testthat); library(scRNAseq); source('test-he-organ-atlas.R')

test_that("Ernst spermatogenesis getter works as expected", {
    sce <- ErnstSpermatogenesisData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- ErnstSpermatogenesisData(method="Cellranger")
    expect_true(ncol(sce)!=ncol(sce2))
})
