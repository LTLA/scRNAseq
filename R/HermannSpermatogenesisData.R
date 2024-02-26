#' Obtain the Hermann spermatogenesis data
#'
#' Obtain the mouse spermatogenesis single-cell RNA-seq data from Hermann et al. (2018).
#'
#' @param strip Logical scalar indicating whether to strip the \code{.X} notation from the row names.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' Only used if \code{strip=TRUE}.
#' @param legacy Logical scalar indicating whether to pull data from ExperimentHub.
#' By default, we use data from the gypsum backend.
#'
#' @details
#' Column metadata contains cell types provided by the data generators
#' at \url{https://data.mendeley.com/datasets/kxd5f8vpt4/1#file-fe79c10b-c42e-472e-9c7e-9a9873d9b3d8}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/hermann-spermatogenesis}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with two matrices, containing spliced and unspliced counts, respectively.
#'
#' @author Charlotte Soneson
#'
#' @references
#' Hermann B.P. et al. (2018).
#' The Mammalian Spermatogenesis Single-Cell Transcriptome, from Spermatogonial Stem Cells to Spermatids. 
#' \emph{Cell Rep.} 25(6), 1650-1667.e8. 
#'
#' @examples
#' sce <- HermannSpermatogenesisData()
#' 
#' @export
HermannSpermatogenesisData <- function(strip=FALSE, location=TRUE, legacy=FALSE) {
    if (!legacy) {
        sce <- fetchDataset("splicing-demonstration-2020", "2023-12-20", realize.assays=TRUE)
    } else {
        version <- "2.4.0"
        sce <- .create_sce(file.path("hermann-spermatogenesis", version), 
            assays = c("spliced", "unspliced"), has.rowdata=FALSE)
    }

    if (strip) {
        rownames(sce) <- sub("\\..*", "", rownames(sce))

        # We know we mapped it to Gencode vM24 --> Ensembl 99, 
        # hence the hard-coded AnnotationHub ID.
        sce <- .define_location_from_ensembl(sce, ahub.id="AH78811", location=location)
    }

    sce 
}
