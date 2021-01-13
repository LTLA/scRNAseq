# library(testthat); library(scRNAseq); source('test-bacher-tcell.R')

test_that("Bacher T cell getter works as expected", {
    sce <- BacherTCellData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- BacherTCellData(filtered=FALSE)
    expect_true(ncol(sce) < ncol(sce2))

    sce <- BacherTCellData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
