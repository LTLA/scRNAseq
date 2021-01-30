# library(testthat); library(scRNAseq); source('test-zhao-immune-liver.R')

test_that("ZhaoImmuneLiverData works as expected", {
    sce <- ZhaoImmuneLiverData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- ZhaoImmuneLiverData(filter=TRUE)
    expect_s4_class(ncol(sce) > ncol(sce2))
    expect_null(sce2$retained)
})
