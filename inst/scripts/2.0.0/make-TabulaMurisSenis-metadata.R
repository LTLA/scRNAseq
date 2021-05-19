require(tidyverse)

droplet_tiss <- c("Bladder",   
                  "Fat",                                     
                  "Heart_and_Aorta",                             
                  "Kidney",                           
                  "Large_Intestine",                      
                  "Limb_Muscle",                        
                  "Liver",                   
                  "Lung",                               
                  "Mammary_Gland",                                
                  "Marrow",                                   
                  "Pancreas",                             
                  "Skin",                        
                  "Spleen",            
                  "Thymus",       
                  "Tongue",       
                  "Trachea")
facs_tiss    <- c("Aorta",                     
                  "BAT",                                  
                  "Bladder",                                  
                  "Brain_Myeloid",                 
                  "Brain_Non-Myeloid",                 
                  "Diaphragm",                     
                  "GAT",                 
                  "Heart",                           
                  "Kidney",                                 
                  "Large_Intestine",                       
                  "Limb_Muscle",                    
                  "Liver",                    
                  "Lung",              
                  "Mammary_Gland",     
                  "Marrow",              
                  "MAT",               
                  "Pancreas",          
                  "SCAT",              
                  "Skin",              
                  "Spleen",            
                  "Tongue",            
                  "Trachea")

datasets     <- c(droplet_tiss, facs_tiss)
technologies <- list()
technologies$droplet<- droplet_tiss
technologies$facs<-facs_tiss
data_types    <- c('counts', 'coldata', 'rowdata')
metadata<- NULL

for (tech_type in names(technologies)){
    for (tissue in technologies[[tech_type]]) {
        for (data in data_types) {
        
          if (tech_type == 'droplet'){descr <- 'Droplet-based (10x Genomics)'} else{descr <- 'FACS sorting - Smart-seq2 technology'}
          h5ad_file_name <- sprintf('tabula-muris-senis-%s-processed-official-annotations-%s%s', tech_type, tissue,'.h5ad')
    
    
          meta <- data.frame(Title= sprintf("Tabula-Muris-Senis-%s-%s-%s", tech_type, tissue, data),
                             Description=sprintf('%s single-cell RNA-seq data from Tabula Muris Senis project 2020', descr),
                             RDataPath= file.path("scRNAseq", sprintf("Tabula-Muris-Senis-%s",tech_type), "2.0.0"),
                             BiocVersion="3.12",
                             Genome = 'GRCm38',
                             Species = 'Mus musculus',
                             TaxonomyId = '10090',
                             SourceType        = 'HDF5',
                             SourceUrl         = 'https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728',
                             SourceVersion     = h5ad_file_name,
                             DataProvider      = 'FigShare',
                             SourceLastModifiedDate = as.POSIXct('2020-07-15 20:17 CEST'),
                             DispatchClass = 'Rds',
                             RDataDateAdded = as.POSIXct(Sys.time(), "CET"),
                             RDataClass=if (grepl("counts", target)) "dgCMatrix" else "DFrame",
                             Maintainer = 'Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>'
                             )
  
    if (is.null(metadata)){metadata <- meta} else {metadata <- rbind(metadata, meta)}
  }}}

  
# Write the data out and put in the inst/extdata directory.
write.csv(metadata, file="metadata-TabulaMurisSenis.csv", row.names=FALSE)