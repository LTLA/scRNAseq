# library(testthat); library(scRNAseq); source('test-zeisel-nervous.R')

test_that("Zeisel nervous getter works as expected", {
    sce <- ZeiselNervousData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- ZeiselNervousData(location=FALSE)
    expect_true(all(lengths(sce)==0))
})
