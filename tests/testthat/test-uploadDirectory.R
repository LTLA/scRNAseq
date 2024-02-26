# This checks that uploadDirectory works as expected.
# library(testthat); library(scRNAseq); source("test-uploadDirectory.R")

sce <- SingleCellExperiment(list(counts=matrix(rpois(1000, lambda=1), 100, 10)))
sce$foo <- sample(letters, 10)
rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
colnames(sce) <- head(LETTERS, 10)

test_that("file listing works as expected without any special elements", {
    meta <- list(
        title="My dataset",
        description="This is my dataset",
        taxonomy_id="10090",
        genome="GRCh38",
        sources=list(list(provider="GEO", id="GSE12345")),
        maintainer_name="Shizuka Mogami",
        maintainer_email="mogami.shizuka@765pro.com"
    )

    tmp <- tempfile()
    saveDataset(sce, tmp, meta)

    listing <- scRNAseq:::list_files(tmp)
    expect_identical(nrow(listing$links), 0L)
    expect_identical(sort(listing$files), sort(list.files(tmp, recursive=TRUE)))
})

test_that("file listing works with ReloadedArrays", {
    cache <- tempfile()
    tmp0 <- file.path(cache, gypsum:::BUCKET_CACHE_NAME, "scRNAseq", "test", "foo")
    dir.create(dirname(tmp0), recursive=TRUE)
    alabaster.base::saveObject(sce, tmp0)

    meta <- list(
        title="My dataset",
        description="This is my dataset",
        taxonomy_id=c("10090", "9606"),
        genome=c("GRCh38", "GRCm38"),
        sources=list(list(provider="GEO", id="GSE12345"), list(provider="DOI", id="123asd/231.123")),
        maintainer_name="Kaori Sakuramori",
        maintainer_email="sakuramori.kaori@765pro.com"
    )

    sce2 <- sce
    assay(sce2, "foobar", withDimnames=FALSE) <- alabaster.matrix::ReloadedArray(path=file.path(tmp0, "assays", "0"), seed=matrix(0, 100, 10))

    tmp <- tempfile()
    saveDataset(sce2, tmp, meta)
    expect_error(scRNAseq:::list_files(tmp), "failed to convert")

    listing <- scRNAseq:::list_files(tmp, cache=cache)
    expect_identical(sort(listing$links$to.path), sort(c("assays/0/array.h5", "assays/0/OBJECT")))
    expect_identical(sort(listing$links$from.path), sort(c("assays/1/array.h5", "assays/1/OBJECT")))
    expect_identical(sort(c(listing$files, listing$links$from.path)), sort(list.files(tmp, recursive=TRUE)))
})

test_that("the actual upload works correctly", {
    gh_token <- gypsum::accessToken(request=NULL)
    skip_if(is.null(gh_token))

    meta <- list(
        title="My dataset",
        description="This is my dataset",
        taxonomy_id="10090",
        genome="GRCh38",
        sources=list(list(provider="GEO", id="GSE12345")),
        maintainer_name="Shizuka Mogami",
        maintainer_email="mogami.shizuka@765pro.com"
    )

    tmp <- tempfile()
    saveDataset(sce, tmp, meta)

    version <- as.character(Sys.Date())
    uploadDirectory(tmp, "test", version, probation=TRUE)
    on.exit(gypsum::rejectProbation("scRNAseq", "test", version))

    cache <- tempfile() # use a different cache to avoid pulling down the test.
    roundtrip <- fetchDataset("test", version, cache=cache)
    expect_identical(colData(roundtrip), colData(sce))
    expect_identical(as.matrix(counts(roundtrip)), assay(sce))
    expect_s4_class(counts(roundtrip, withDimnames=FALSE), "ReloadedMatrix")
})
