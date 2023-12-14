# This tests the dataset fetching.
# library(testthat); library(scRNAseq); source("test-fetchDataset.R")

test_that("fetchDataset works as expected", {
    sce <- fetchDataset("zeisel-brain-2015", "2.17.1")
    expect_s4_class(sce, "SingleCellExperiment")
    expect_s4_class(assay(sce, withDimnames=FALSE), "ScrnaseqMatrix")
})

test_that("fetchMetadata works as expected", {
    meta <- fetchMetadata("zeisel-brain-2015", "2.17.1")
    expect_match(meta$title, "Brain structure")
    expect_identical(meta$taxonomy_id[[1]], "10090")
})
