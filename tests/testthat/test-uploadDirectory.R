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
    expect_identical(nrow(data.frame(listing$links)), 0L)
    files <- data.frame(listing$files)
    expect_identical(sort(files$path), sort(list.files(tmp, recursive=TRUE)))
})

test_that("file listing works with ScrnaseqMatrices", {
    tmp0 <- tempfile()
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
    assay(sce2, "foobar", withDimnames=FALSE) <- ScrnaseqArray(
        name="test", 
        version="foo", 
        path="bar", 
        cached=file.path(tmp0, "assays", "0"),
        seed=matrix(0, 100, 10)
    )

    tmp <- tempfile()
    saveDataset(sce2, tmp, meta)
    listing <- scRNAseq:::list_files(tmp)
    links <- data.frame(listing$links)
    expect_identical(links$to.path, c("bar/array.h5", "bar/OBJECT"))
    expect_identical(links$from.path, c("assays/1/array.h5", "assays/1/OBJECT"))

    files <- data.frame(listing$files)
    expect_identical(sort(c(files$path, links$from.path, "assays/1/_link")), sort(list.files(tmp, recursive=TRUE)))

    # Works in the pathological case where path = NULL. This shouldn't happen
    # as there should always be a path to an SE's assay, but you never know.
    {
        sce2 <- sce
        assay(sce2, "foobar", withDimnames=FALSE) <- ScrnaseqArray(
            name="test", 
            version="foo", 
            path=NULL, 
            cached=file.path(tmp0, "assays", "0"),
            seed=matrix(0, 100, 10)
        )

        tmp <- tempfile()
        saveDataset(sce2, tmp, meta)
        listing <- scRNAseq:::list_files(tmp)
        links <- data.frame(listing$links)
        expect_identical(links$to.path, c("array.h5", "OBJECT"))
        expect_identical(links$from.path, c("assays/1/array.h5", "assays/1/OBJECT"))
    }
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

    roundtrip <- fetchDataset("test", version)
    expect_identical(colData(roundtrip), colData(sce))
    expect_identical(as.matrix(counts(roundtrip)), assay(sce))
    expect_s4_class(counts(roundtrip, withDimnames=FALSE), "ScrnaseqMatrix")
})
