# library(testthat); library(scRNAseq); source("test-fletcher-olfactory.R")

test_that("FletcherOlfactoryData works as expected", {
    sce <- FletcherOlfactoryData(ensembl=TRUE)
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- FletcherOlfactoryData(filtered=FALSE)
    expect_true(ncol(sce2) > ncol(sce))
    expect_true(nrow(sce2) > nrow(sce))
})
