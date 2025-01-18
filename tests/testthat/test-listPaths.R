# library(scRNAseq); library(testthat); source("test-listPaths.R")

test_that("listPaths works correctly with names", {
    he.listing <- listPaths("he-organs-2020", "2023-12-21")
    expect_gt(length(he.listing), 2)
    expect_true("blood" %in% he.listing)

    z.listing <- listPaths("zeisel-brain-2015", "2023-12-14") # no subdatasets
    expect_identical(length(z.listing), 1L)
    expect_true(is.na(z.listing))
})

test_that("listPaths works correctly with metadata", {
    he.listing <- listPaths("he-organs-2020", "2023-12-21", include.metadata=TRUE)
    expect_gt(nrow(he.listing), 2)
    expect_true("bladder" %in% he.listing$path)

    z.listing <- listPaths("zeisel-brain-2015", "2023-12-14", include.metadata=TRUE) 
    expect_identical(nrow(z.listing), 1L)
    expect_true(is.na(z.listing$path))
})
