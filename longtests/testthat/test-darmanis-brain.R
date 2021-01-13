# library(testthat); library(scRNAseq); source('test-darmanis-brain.R')

test_that("Darmanis brain getter works as expected", {
    sce <- DarmanisBrainData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce <- DarmanisBrainData(ensembl=TRUE)
    expect_true(all(grepl("^ENSG", rownames(sce))))
})
