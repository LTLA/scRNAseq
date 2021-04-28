# library(testthat); library(scRNAseq); source('test-bunis-hspc.R')

test_that("Bunis HSPC getter works as expected", {
  sce <- BunisHSPCData()
  expect_s4_class(sce, "SingleCellExperiment")
  
  sce2 <- BunisHSPCData(filtered = "cells")
  expect_s4_class(sce2, "SingleCellExperiment")
  
  sce3 <- BunisHSPCData(filtered=FALSE)
  
  # Checks cell filtering and that colData is added all at once
  expect_true( nrow(colData(sce)) < nrow(colData(sce2)) )
  expect_true( nrow(colData(sce2)) < nrow(colData(sce3)) )
  
  expect_true(all(grepl("^ENSG", rownames(sce))))
})
