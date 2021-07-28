# library(testthat); library(scRNAseq); source('test-ledergor-myeloma.R')

test_that("Ledegor myeloma getter works as expected", {
    sce <- LedergorMyelomaData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- LedergorMyelomaData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
