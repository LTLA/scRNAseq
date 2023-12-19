# This checks that saveDataset works as expected.
# library(testthat); library(scRNAseq); source("test-polishDataset.R")

mat <- matrix(rpois(1000, lambda=1), ncol=10)

test_that("polishDataset strips assay dimnames", {
    rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
    colnames(mat) <- head(LETTERS, 10)
    sce <- SingleCellExperiment(list(counts=mat))

    y <- polishDataset(sce)
    expect_identical(dimnames(y), dimnames(mat))
    expect_null(dimnames(assay(y, withDimnames=FALSE)))
})

test_that("polishDataset strips reduced dimension names", {
    sce <- SingleCellExperiment(list(counts=mat))
    rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(mat)))
    colnames(sce) <- head(LETTERS, 10)

    pca <- matrix(rnorm(50), ncol=5)
    colnames(pca) <- sprintf("PC%i", seq_len(5))
    rownames(pca) <- seq_len(nrow(pca))
    reducedDim(sce, "PCA", withDimnames=FALSE) <- pca

    y <- polishDataset(sce)
    rdim <- reducedDim(y, "PCA", withDimnames=FALSE)
    expect_null(rownames(rdim))
    expect_identical(colnames(rdim), colnames(pca))
})

test_that("polishDataset strips reduced dimension names", {
    sce <- SingleCellExperiment(list(counts=mat))
    rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(mat)))
    colnames(sce) <- head(LETTERS, 10)

    sub <- sce[1:10,]
    rownames(sub) <- paste0("X-", rownames(sub))
    altExp(sce, "ERCC", withDimnames=FALSE) <- sub

    y <- polishDataset(sce)
    ae <- altExp(y, "ERCC", withDimnames=FALSE)
    expect_identical(rownames(ae), rownames(sub))
    expect_null(colnames(ae))
})

library(DelayedArray)
test_that("polishDataset converts dense to sparse", {
    mat <- matrix(rpois(1000, lambda=0.2), ncol=10)
    sce <- SingleCellExperiment(list(counts=mat))

    y <- polishDataset(sce)
    expect_true(is_sparse(assay(y)))

    # No conversion.
    y <- polishDataset(sce, reformat.assay.by.density=0)
    expect_false(is_sparse(assay(y)))
    y <- polishDataset(sce, reformat.assay.by.density=NULL)
    expect_false(is_sparse(assay(y)))
})

library(Matrix)
test_that("polishDataset converts sparse to dense", {
    mat <- matrix(rpois(1000, lambda=3), ncol=10)
    sce <- SingleCellExperiment(list(counts=as(mat, "dgCMatrix")))

    y <- polishDataset(sce)
    expect_false(is_sparse(assay(y)))

    # No conversion.
    y <- polishDataset(sce, reformat.assay.by.density=1)
    expect_true(is_sparse(assay(y)))
    y <- polishDataset(sce, reformat.assay.by.density=NULL)
    expect_true(is_sparse(assay(y)))
})

test_that("polishDataset attempts integer conversions", {
    mat <- matrix(rpois(1000, lambda=3), ncol=10) * 1.0
    sce <- SingleCellExperiment(list(counts=mat))
    expect_identical(type(assay(sce)), "integer")

    y <- polishDataset(sce)
    expect_identical(type(assay(y)), "integer")

    # Handles the conversion safely.
    mat <- matrix(rpois(1000, lambda=0.1), ncol=10) 
    as(mat, "dgCMatrix")))

})

