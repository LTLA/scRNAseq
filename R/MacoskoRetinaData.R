#' Obtain the Macosko retina data
#'
#' Obtain the mouse retina single-cell RNA-seq data from Macosko et al. (2016).
#'
#' @param ensembl Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#'
#' @details
#' Column metadata contains the cluster identity as reported in the paper.
#' Note that some cells will have \code{NA} identities as they are present in the count matrix but not in the metadata file.
#' These are presumably low-quality cells that were discarded prior to clustering.
#'
#' If \code{ensembl=TRUE}, the gene symbols are converted to Ensembl IDs in the row names of the output object.
#' Rows with missing Ensembl IDs are discarded, and only the first occurrence of duplicated IDs is retained.
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#' Note that this is only performed if \code{ensembl=TRUE}.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/macosko-retina}.
#'
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' Macosko E et al. (2016). 
#' Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. 
#' \emph{Cell} 161(5), 1202-1214.
#'
#' @examples
#' sce <- MacoskoRetinaData()
#' 
#' @export
MacoskoRetinaData <- function(ensembl=FALSE, location=TRUE) {
    version <- "2.0.0"
    sce <- .create_sce(file.path("macosko-retina", version), has.rowdata=FALSE)

    if (ensembl) {
        # For some bizarre reason, this dataset has all-caps symbols... for mice.
        # So we need to do some custom work to ensure that this converts properly.
        tag <- "AH73905"
        edb <- AnnotationHub::AnnotationHub()[[tag]]

        anno <- AnnotationDbi::select(edb, keys=AnnotationDbi::keys(edb), keytype="GENEID", columns="SYMBOL")
        all.symbols <- tolower(anno$SYMBOL)
        cur.symbols <- tolower(rownames(sce)) 

        ensid <- anno$GENEID[match(cur.symbols, all.symbols)]
        keep <- !is.na(ensid) & !duplicated(ensid)

        sce <- sce[keep, ]
        rownames(sce) <- ensid[keep]

        sce <- .define_location_from_ensembl(sce, species="Mm", location=location)
    }

    sce
}
