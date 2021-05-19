#' Obtain the Tabula Muris Senis droplet data
#'
#' Obtain droplet-based (10x Genomics) single-cell RNA-seq data from Tabula Muris Senis project 2020.
#' This is a comprehensive analysis of the aging dynamics across the mouse lifespan.
#'
#' @param tissue   A string specifying which datasets should be obtained.
#' @param age      A character vector specifying age samples to retrieve given a specific tissue.
#' @param ensembl  Logical scalar indicating whether the output row names should contain Ensembl identifiers.
#' @param location Logical scalar indicating whether genomic coordinates should be returned.
#' 
#' @details
#' If age in \code{age} is not specified the all ages  for that tissue will be retrieved.
#' Cells belonging to different age are recognizable from column 'age' in metadata.
#' Be aware that not all tissues have samples for all ages ("1m", "3m", "18m",  "21m", "24m","30m").
#'
#' If \code{location=TRUE}, the coordinates of the Ensembl gene models are stored in the \code{\link{rowRanges}} of the output.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#' Specific resources can be retrieved by searching for \code{scRNAseq/TabulaMurisSenisDroplet}.
#' 
#' @return A \linkS4class{SingleCellExperiment} object with a single matrix of UMI counts.
#'
#' @author Aaron Lun
#'
#' @references
#' The Tabula Muris Consortium (2020).
#' A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. 
#' \emph{Nature} 583, 590–595(2020)
#'
#' @examples
#' sce <- TabulaMurisSenisDropletData( tissue=c("Bladder"),
#'                                     age=c("1m", "3m", "18m",  "24m"),
#'                                     ensembl=FALSE,
#'                                     location=TRUE)
#' 
#' @export
#' @importFrom SummarizedExperiment rowData<- rowData
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom tidyverse tribble

TabulaMurisSenisDropletData <- function( tissue=c("Bladder"),
                                         age=c("1m", "3m", "18m",  "24m"),
                                         ensembl=FALSE,
                                         location=TRUE)
{


                          
  Overview <- tribble(
    ~"Tissue",              ~'1m',  ~'3m',  ~'18m',  ~'21m',  ~'24m',  ~'30m',
    #--------------------|--------|--------|--------|--------|--------|--------|
    "Bladder",              "X",      "X",    "X",     "...",   "X",     "...",     
    "Fat",                  "...",    "...",  "X",     "X",     "...",   "X",                    
    "Heart_and_Aorta",      "X",      "X",    "X",     "X",     "X",     "X",                        
    "Kidney",               "X",      "X",    "X",     "X",     "X",     "X",             
    "Large_Intestine",      "...",    "...",  "...",   "...",   "...",   "X",                 
    "Limb_Muscle",          "X",      "X",    "X",     "X",     "X",     "X",               
    "Liver",                "X",      "X",    "X",     "X",     "X",     "X",    
    "Lung",                 "X",      "X",    "X",     "X",     "...",   "X",               
    "Mammary_Gland",        "...",    "X",    "X",     "X",     "...",   "...",                            
    "Marrow",               "X",      "X",    "X",     "X",     "X",     "X",                     
    "Pancreas",             "...",    "...",  "X",     "X",     "X",     "X",                 
    "Skin",                 "...",    "...",  "X",     "X",     "...",   "...",          
    "Spleen",               "X",      "X",    "X",     "X",     "X",     "X",          
    "Thymus",               "...",    "X",    "X",     "X",     "X",     "...",       
    "Tongue",               "X",      "X",    "X",     "...",   "X",     "...",        
    "Trachea",              "...",    "X",    "...",   "...",   "...",   "...",       
     )                      

  if (length(tissue) > 1) {print('Select only one tissue from the following list:')
    print(Overview)}
  else{
                           
  version <- "2.0.0"
  host    <- file.path("Tabula-Muris-Senis-droplet", version)

  sce     <- .create_sce(host, suffix=tissue)
  
  # Select the age of interest
  if (!is_null(age)){
    age <- age[age %in% levels(colData(sce)$age)]
    sce <- sce[,colData(sce)$age %in% age]
  }

  .convert_to_ensembl(sce, 
                      symbols=rownames(sce), 
                      species="Mm",
                      ensembl=ensembl,
                      location=location)
  }
}