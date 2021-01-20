# library(testthat); library(scRNAseq); source('test-bhaduri-organoid.R')

test_that("Bhaduri organoid getter works as expected", {
    sce <- BhaduriOrganoidData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- BhaduriOrganoidData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
