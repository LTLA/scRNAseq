library(scRNAseq)
df = enhanceListDatasets()
library(hwriter)
hwrite(df, "new.html")
z = readLines("new.html")
z = z[-c(1:13)] # remove header material
z = c("# scRNAseq package", 
#This package provides functions for retrieving SingleCellExperiment",
representations of data from 55 single-cell RNA-seq experiments.",
"",
"Functions to call after running `library(scRNAseq)` are",
"listed in the 'call' column of the table:", z)
writeLines(z, "README.md")

