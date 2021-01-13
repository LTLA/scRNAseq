# library(testthat); library(scRNAseq); source('test-zhong-prefrontal.R')

test_that("Zhong prefrontal getter works as expected", {
    sce <- ZhongPrefrontalData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- ZhongPrefrontalData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
