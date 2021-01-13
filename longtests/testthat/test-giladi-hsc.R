# library(testthat); library(scRNAseq); source('test-giladi-hsc.R')

test_that("Giladi HSC getter works as expected", {
    sce <- GiladiHSCData()
    expect_s4_class(sce, "SingleCellExperiment")

    sce2 <- GiladiHSCData(mode="rna")
    expect_true(ncol(sce2) > ncol(sce))

    sce.filt <- GiladiHSCData(mode='rna', filtered=TRUE)
    expect_identical(colnames(sce.filt), colnames(sce2)[sce2$retained])

    sce3 <- GiladiHSCData(mode="crispr")
    expect_true(ncol(sce3) > ncol(sce))

    sceE <- GiladiHSCData(mode="rna", ensembl=TRUE)
    expect_true(all(grepl("^ENSMUSG", rownames(sceE))))
})
