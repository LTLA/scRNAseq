#' Obtain the Stoeckius cell hashing data
#'
#' Obtain the (mostly human) cell hashing single-cell RNA-seq data from Stoeckius et al. (2018).
#'
#' @param type String specifying the dataset to obtain.
#' @param mode String specifying the data modalities to obtain, see Details.
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' @param strip.metrics Logical scalar indicating whether quality control metrics should be removed from the HTO/ADT counts.
#'
#' @details
#' When \code{type="pbmc"}, the \code{mode} can be one or more of:
#' \itemize{
#' \item \code{"human"}, the RNA counts for human genes.
#' \item \code{"mouse"}, the RNA counts for mouse genes.
#' Present as the PBMC dataset is actually a mixture of human PBMCs and unlabelled mouse cells.
#' \item \code{"hto"}, the HTO counts.
#' \item \code{"adt1"}, counts for the first set of ADTs (immunoglobulin controls).
#' \item \code{"adt2"}, counts for the second set of ADTs (cell type-specific markers).
#' }
#'
#' When \code{type="mixed"}, the \code{mode} can be one or more of:
#' \itemize{
#' \item \code{"rna"}, the RNA counts for the genes;
#' \item \code{"hto"}, the HTO counts.
#' }
#'
#' If \code{ensembl=TRUE}, gene symbols for the RNA counts are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE} and only for the RNA counts.
#'
#' For the HTO and ADT matrices, some rows correspond to quality control metrics.
#' If \code{strip.metrics=TRUE}, these rows are removed so that only data for actual HTOs or ADTs are present.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/nestorowa-hsc}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a matrix of UMI counts corresponding to the first \code{mode},
#' plus any number of alternative Experiments containing the remaining \code{mode}s.
#' If multiple \code{mode}s are specified, the output object only contains the intersection of their column names.
#'
#' @author Aaron Lun
#'
#' @references
#' Stoeckius et al. (2018). 
#' Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics.
#' \emph{Genome Biol.} 19, 224.
#'
#' @examples
#' sce <- StoeckiusHashingData()
#' 
#' @export
#' @importFrom SingleCellExperiment altExps<- 
StoeckiusHashingData <- function(type=c("pbmc", "mixed"), mode=NULL,
    ensembl=FALSE, location=TRUE, strip.metrics=TRUE) 
{
    version <- "2.4.0"

    type <- match.arg(type)
    if (is.null(mode)) {
        if (type=="pbmc") {
            mode <- c("human", "mouse", "hto")
        } else {
            mode <- c("rna", "hto")
        }
    }

    if (type=="pbmc") {
        acceptable <- c("human", "mouse", "hto", "adt1", "adt2")
    } else {
        acceptable <- c("rna", "hto")
    }
    mode <- match.arg(mode, acceptable, several.ok=TRUE)

    collated <- list()
    for (m in mode) {
        collated[[m]] <- .create_sce(file.path("stoeckius-hashing", version), 
            suffix=paste0(type, "-", m), has.coldata=FALSE)
    }

    if (length(collated) > 1) {
        common <- Reduce(intersect, lapply(collated, colnames))
        for (m in seq_along(collated)) {
            collated[[m]] <- collated[[m]][,common]
        }
    }

    # Converting everything to Ensembl and adding locations.
    if (type=="pbmc") {
        convertible <- c(human="Hs", mouse="Mm")
    } else {
        convertible <- c(rna="Hs")
    }

    for (i in intersect(names(convertible), names(collated))) {
        collated[[i]] <- .convert_to_ensembl(collated[[i]],
            species=converted[i],
            symbols=rownames(collated[[i]]),
            ensembl=ensembl,
            location=location)
    }

    # Stripping metrics from the ADT matrices.
    if (strip.metrics) {
        strippable <- c("hto", "adt1", "adt2")
        metrics <- c("no_match", "ambiguous", "total_reads", "bad_struct")
        for (i in intersect(names(collated), strippable)) {
            collated[[i]] <- collated[[i]][setdiff(rownames(collated[[i]]), metrics),]
        }
    }

    primary <- collated[[1]]
    altExps(primary) <- collated[-1]
    primary
}
