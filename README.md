# Single-cell RNA-seq datasets for Bioconductor

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html)|[![Release OK](https://bioconductor.org/shields/build/release/data-experiment/scRNAseq.svg)](http://bioconductor.org/checkResults/release/data-experiment-LATEST/scRNAseq/)|
|[BioC-devel](https://bioconductor.org/packages/devel/data/experiment/html/scRNAseq.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/data-experiment/scRNAseq.svg)](http://bioconductor.org/checkResults/devel/data-experiment-LATEST/scRNAseq/)|

This package provides single-cell (mostly RNA-seq) datasets for convenient use by other Bioconductor packages and workflows.
Each dataset is loaded as a [`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment) that is immediately ready for further analysis.
To get started, install the package and its dependencies from Bioconductor:

```r
# install.packages("BiocManager")
BiocManager::install("scRNAseq")
```

Find datasets of interest:

```r
surveyDatasets()
searchDatasets("brain")
```

Fetch a dataset as a `SingleCellExperiment`:

```r
fetchDataset("zeisel-brain-2015", version="2023-12-14")
fetchDataset("baron-pancreas-2016", "2023-12-14", path="human")
```

And add your own datasets to enable re-use by the wider Bioconductor community.

Check out the [user's guide](https://bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html) for more details.

## Maintainer notes

If someone wants to contribute a new dataset, they should follow the instructions in the user's guide.
This requires a bit of effort on our (i.e., the package maintainers') part to process a new submission.

Once a prospective uploader has prepared their `SingleCellExperiment` and its associated metadata,
they can be given temporary upload permissions for, e.g., a week, by calling:

```r
gypsum::setPermissions("scRNAseq", uploaders=list(
    list(
        id="GITHUB_LOGIN", 
        asset="NAME_OF_THE_DATASET_THEY_WANT_TO_UPLOAD",
        version="VERSION_THEY_WANT_TO_UPLOAD",
        until=Sys.time() + 7 * 24 * 3600
    )
)
```

Once the upload is complete, we pull down the dataset for review.

```r
cache <- tempfile()
dest <- gypsum::saveVersion(
    "scRNAseq", 
    asset="NAME_OF_THE_DATASET_THEY_WANT_TO_UPLOAD",
    version="VERSION_THEY_WANT_TO_UPLOAD",
    cache=cache
)

# Check that the saved object is valid. You might need to do this for each
# subdirectory if they saved multiple objects in a single dataset.
alabaster.base::validateObject(dest)

# You can also just try loading it for inspection in the R session.
alabaster.base::readObject(dest)

# Check that the metadata is valid.
lines <- readLines(file.path(dest, "_bioconductor.json"))
gypsum::validateMetadata(paste(lines, collapse="\n"))
```

If everything looks okay, we can approve the probational dataset.
Otherwise we reject it.

```r
# Okay.
gypsum::approveProbation(
    "scRNAseq", 
    asset="NAME_OF_THE_DATASET_THEY_WANT_TO_UPLOAD",
    version="VERSION_THEY_WANT_TO_UPLOAD"
)

# Not okay.
gypsum::rejectProbation(
    "scRNAseq", 
    asset="NAME_OF_THE_DATASET_THEY_WANT_TO_UPLOAD",
    version="VERSION_THEY_WANT_TO_UPLOAD"
)
```
