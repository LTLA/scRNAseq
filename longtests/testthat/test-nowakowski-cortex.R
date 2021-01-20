# library(testthat); library(scRNAseq); source('test-nowakowski-cortex.R')

test_that("Nowakowski cortex getter works as expected", {
    sce <- NowakowskiCortexData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- NowakowskiCortexData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
