#' Polish dataset for saving
#'
#' Prepare a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} to be saved with \code{\link{saveDataset}}.
#' This performs some minor changes to improve storage efficiency.
#'
#' @param x A \linkS4class{SummarizedExperiment} or one of its subclasses.
#' @param strip.names Logical scalar indicating whether to strip redundant names from internal objects, 
#' e.g., dimnames of the assays, row names of reduced dimensions, column names of alternative experiments.
#' This saves some space in the on-disk representation.
#' @param reformat.assay.by.density Numeric scalar indicating whether to optimize assay formats based on the density of non-zero values.
#' Assays with densities above this number are converted to ordinary dense arrays (if they are not already), while those with lower densities are converted to sparse matrices.
#' This can be disabled by setting it to \code{NULL}.
#' @param attempt.integer.conversion Logical scalar indicating whether to convert double-precision assays containing integer values to actually have the integer type.
#' This can improve efficiency of downstream applications by avoiding the need to operate in double precision.
#' @param forbid.nested.altexp Logical scalar indicating whether nested alternative experiments (i.e., alternative experiments of alternative experiments) should be forbidden.
#' This defaults to \code{TRUE} as nested alternative experiments are usually the result of some mistake in \code{\link{altExp}} preparation.
#' @param remove.altexp.coldata Logical scalar indicating whether column data for alternative experiments should be removed.
#' This defaults to \code{TRUE} as the alternative experiment column data is usually redundant with that of the main experiment.
#'
#' @return A modified copy of \code{x}.
#'
#' @author Aaron Lun
#'
#' @examples
#' mat <- matrix(rpois(1000, lambda=0.2), 100, 10) * 1.0
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#' colnames(mat) <- head(LETTERS, 10)
#'
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(list(counts=mat))
#' str(assay(sce, withDimnames=FALSE))
#'
#' polished <- polishDataset(sce)
#' str(assay(polished, withDimnames=FALSE))
#'
#' @export
polishDataset <- function(x, strip.inner.names=TRUE, reformat.assay.by.density=0.3, attempt.integer.conversion=TRUE, remove.altexp.coldata=TRUE, forbid.nested.altexp=TRUE) {
    .polish_dataset(x, 
        strip.inner.names=strip.inner.names, 
        reformat.assay.by.density=reformat.assay.by.density, 
        attempt.integer.conversion=attempt.integer.conversion, 
        remove.altexp.coldata=remove.altexp.coldata, 
        forbid.nested.altexp=forbid.nested.altexp
    )
}

#' @importFrom BiocGenerics mean
#' @importFrom DelayedArray is_sparse type<- DelayedArray
#' @importClassesFrom SparseArray SVT_SparseMatrix
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment reducedDim reducedDimNames reducedDim<- altExp altExpNames altExp<-
.polish_dataset <- function(x, strip.inner.names, reformat.assay.by.density, attempt.integer.conversion, remove.altexp.coldata, forbid.nested.altexp, level=0) {
    for (i in assayNames(x)) {
        ass <- assay(x, i, withDimnames=FALSE)

        if (!is.null(reformat.assay.by.density)) {
            density <- mean(ass != 0 | is.na(ass)) # NAs are included as non-zero.
            if (density < reformat.assay.by.density) {
                if (!is_sparse(ass)) {
                    ass <- as(ass, "SVT_SparseMatrix")
                }
            } else {
                if (is_sparse(ass)) {
                    ass <- as.array(ass)
                }
            }
        }

        if (attempt.integer.conversion) {
            if (type(ass) == "double" && !any(ass %% 1 != 0, na.rm=TRUE)) {
                converter <- selectMethod(`type<-`, class(ass), optional=TRUE)
                if (is.null(converter)) {
                    ass <- DelayedArray(ass)
                    type(ass) <- "integer"
                } else {
                    ass <- converter(ass, "integer")
                }
            }
        }

        if (strip.inner.names) {
            dimnames(ass) <- NULL
        }

        assay(x, i, withDimnames=FALSE) <- ass
    }

    if (is(x, "SingleCellExperiment")) {
        if (strip.inner.names) {
            for (i in reducedDimNames(x)) {
                red <- reducedDim(x, i, withDimnames=FALSE)
                rownames(red) <- NULL
                reducedDim(x, i, withDimnames=FALSE) <- red
            }
        }

        for (i in altExpNames(x)) {
            if (forbid.nested.altexp && level > 0L) {
                stop("nested alternative experiments are forbidden")
            }

            alt <- altExp(x, i, withDimnames=FALSE)
            if (remove.altexp.coldata) {
                colData(alt) <- colData(alt)[,0]
            }
            if (strip.inner.names) {
                colnames(alt) <- NULL
            }

            alt <- .polish_dataset(alt, 
                strip.inner.names=strip.inner.names,
                reformat.assay.by.density=reformat.assay.by.density, 
                attempt.integer.conversion=attempt.integer.conversion,
                remove.altexp.coldata=remove.altexp.coldata,
                forbid.nested.altexp=forbid.nested.altexp,
                level=level + 1L
            )

            altExp(x, i, withDimnames=FALSE) <- alt
        }
    }

    x
}
