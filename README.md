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
Check out the [user's guide](https://bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html) for more information.
