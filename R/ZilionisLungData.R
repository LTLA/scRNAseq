#' Obtain the Zilionis lung cancer data
#'
#' Obtain the human/mouse lung cancer single-cell RNA-seq data from Zilionis et al. (2019).
#'
#' @param which String specifying the species to get data for.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param filter Logical scalar indicating if the filtered subset should be returned.
#'
#' @details
#' Column metadata is provided and contains information on the library, donor ID/animal ID, replicate and tissue.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#' 
#' If \code{filter=TRUE}, only cells that have been used in the original analysis are returned.
#' In this case, the returned object will also contain coordinates of a SPRING representation in its \code{\link{reducedDims}}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/zilionis-lung}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of read counts.
#'
#' @author Jens Preussner
#'
#' @references
#' Zilionis R et al. (2019). 
#' Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species. 
#' \emph{Immunity} 50(5), 1317-1334.
#'
#' @examples
#' sce.human <- ZilionisLungData()
#'
#' sce.mouse <- ZilionisLungData("mouse")
#' 
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
ZilionisLungData <- function(which=c("human", "mouse"), ensembl=FALSE, location=TRUE, filter=FALSE) {
  version <- "2.4.0"
  which <- match.arg(which)
  sce <- .create_sce(file.path("zilionis-lung", version), has.rowdata=FALSE, suffix=which)
  
  if(filter) {
      # Subset to cells from the original analysis
      keep <- !is.na(colData(sce)$Total.counts)
      sce <- sce[, keep]
      
      # Add spring representation
      which_x <- grep("^x_.*_all|^x$", colnames(colData(sce)))
      which_y <- grep("^y_.*_all|^y$", colnames(colData(sce)))
      spring_rep <- as.matrix(colData(sce)[, c(which_x, which_y)])
      reducedDims(sce, "SPRING") <- spring_rep

      colData(sce) <- colData(sce)[,-c(which_x, which_y)]
  }
  
  .convert_to_ensembl(sce, 
                      ensembl=ensembl,
                      species=if (which=="human") "Hs" else "Mm",
                      symbols=rownames(sce),
                      location=location)
}
