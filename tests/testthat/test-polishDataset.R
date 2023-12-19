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

    # Unless we disable it.
    rownames(mat) <- tolower(rownames(mat))
    colnames(mat) <- tolower(colnames(mat))
    assay(sce, withDimnames=FALSE) <- mat
    y <- polishDataset(sce, strip.inner.names=FALSE)
    expect_identical(dimnames(y), dimnames(sce))
    expect_identical(dimnames(assay(y, withDimnames=FALSE)), dimnames(mat))
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

    # Unless we disable it.
    y <- polishDataset(sce, strip.inner.names=FALSE)
    rdim <- reducedDim(y, "PCA", withDimnames=FALSE)
    expect_identical(dimnames(rdim), dimnames(pca))
})

test_that("polishDataset strips alternative experiment names", {
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

    # Unless we disable it.
    y <- polishDataset(sce, strip.inner.names=FALSE)
    ae <- altExp(y, "ERCC", withDimnames=FALSE)
    expect_identical(rownames(ae), rownames(sub))
    expect_identical(colnames(ae), colnames(sub))
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

    y <- polishDataset(sce)
    a <- assay(y)
    expect_true(is.matrix(a))
    expect_type(a, "integer")

    # Handles the conversion safely.
    mat <- matrix(rpois(1000, lambda=0.1), ncol=10) * 1.0
    sce <- SingleCellExperiment(list(counts=mat))

    y <- polishDataset(sce)
    a <- assay(y)
    expect_s4_class(a, "SVT_SparseMatrix")
    expect_identical(type(a), "integer")

    # Handles the conversion when there's no type<- method:
    mat <- round(as(mat, "dgCMatrix") * 100)
    sce <- SingleCellExperiment(list(counts=mat))

    y <- polishDataset(sce)
    a <- assay(y)
    expect_s4_class(a, "DelayedMatrix")
    expect_identical(type(a), "integer")

    # Doesn't attempt the conversion.
    mat <- matrix(rpois(1000, lambda=3), ncol=10) * 1.5
    sce <- SingleCellExperiment(list(counts=mat))

    y <- polishDataset(sce)
    a <- assay(y)
    expect_true(is.matrix(a))
    expect_type(a, "double")
})

test_that("polishDataset works with NA values in the matrices", {
    mat <- matrix(rpois(1000, lambda=0.1), ncol=10)
    mat[sample(length(mat), 10)] <- NA
    sce <- SingleCellExperiment(list(counts=as(mat, "dgCMatrix")))

    y <- polishDataset(sce)
    expect_true(is_sparse(assay(y)))
    expect_identical(type(assay(y)), "integer")
})

test_that("polishDataset strips out alternative experiment's colData", {
    sce <- SingleCellExperiment(list(counts=mat))
    rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(mat)))
    colnames(sce) <- head(LETTERS, 10)

    sub <- sce[1:10,]
    rownames(sub) <- paste0("X-", rownames(sub))
    sub$foo <- runif(10)
    altExp(sce, "ERCC", withDimnames=FALSE) <- sub

    y <- polishDataset(sce)
    expect_identical(ncol(altExp(y)), 0L)

    # Unless we disable it.
    y <- polishDataset(sce, remove.altexp.coldata=FALSE)
    expect_identical(ncol(altExp(y)), 1L)
})

test_that("polishDataset forbids highly nested altexps", {
    sce0 <- SingleCellExperiment(list(counts=mat))
    sce1 <- SingleCellExperiment(list(counts=mat[1:10,]+1))
    sce2 <- SingleCellExperiment(list(counts=mat[1:5,]+2))
    altExp(sce1, "FOO") <- sce2
    altExp(sce0, "BAR") <- sce1

    expect_error(polishDataset(sce0), "nested")
    y <- polishDataset(sce0, forbid.nested.altexp=FALSE)
    expect_equal(y, sce0)
})
