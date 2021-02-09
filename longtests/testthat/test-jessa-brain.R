# library(testthat); library(scRNAseq); source("test-jessa-brain.R")

test_that("JessaBrainData works as expected", {
    sce <- JessaBrainData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- JessaBrainData(filtered=FALSE)
    expect_true(ncol(sce2) > ncol(sce))
})
