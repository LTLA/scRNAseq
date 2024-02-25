# This tests the dataset fetching.
# library(testthat); library(scRNAseq); source("test-fetchDataset.R")

test_that("fetchDataset works as expected", {
    sce <- fetchDataset("zeisel-brain-2015", "2023-12-14")
    expect_s4_class(sce, "SingleCellExperiment")

    # Correctly creates ReloadedMatrix objects.
    ass <- assay(sce, withDimnames=FALSE)
    expect_s4_class(ass, "ReloadedMatrix")
    expect_true(DelayedArray::is_sparse(ass))
    expect_true(grepl("zeisel-brain-2015", ass@seed@path))
    expect_true(grepl("2023-12-14", ass@seed@path))

    # Works with realization options.
    sce <- fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays=TRUE)
    expect_s4_class(assay(sce, withDimnames=FALSE), "dgCMatrix")
    expect_type(assay(altExp(sce), withDimnames=FALSE), "integer") # also realizes the alternative experiments.
})

test_that("fetchDataset realizes the reduced dimensions", {
    sce <- fetchDataset("aztekin-tail-2019", "2023-12-14", realize.reduced.dims=FALSE)
    expect_s4_class(reducedDim(sce, withDimnames=FALSE), "ReloadedMatrix")

    sce <- fetchDataset("aztekin-tail-2019", "2023-12-14", realize.reduced.dims=TRUE)
    expect_type(reducedDim(sce, withDimnames=FALSE), "double")
})

test_that("fetchMetadata works as expected", {
    meta <- fetchMetadata("zeisel-brain-2015", "2023-12-14")
    expect_match(meta$title, "Brain structure")
    expect_identical(meta$taxonomy_id[[1]], "10090")
})
