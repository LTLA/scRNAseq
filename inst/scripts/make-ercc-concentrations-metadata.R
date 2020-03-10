# mkdir -p scRNAseq/ercc-concentrations/2.2.0/
# cd scRNAseq/ercc-concentrations/2.2.0/
# wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt

write.csv(file = "../extdata/2.2.0/metadata-ercc-concentrations.csv", 
    data.frame(
        Title = "ERCC spike-in concentrations",
        Description = "Concentrations for ERCC spike-in RNA molecules in 2 mixes",
        RDataPath = file.path("scRNAseq", "ercc-concentrations", "2.2.0", "cms_095046.txt"),
        BiocVersion = "3.11",
        Genome = NA,
        SourceType = "TSV",
        SourceUrl = "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/",
        SourceVersion = "cms_095046.txt",
        Species = NA,
        TaxonomyId = NA,
        Coordinate_1_based = NA,
        DataProvider = "Thermo Fisher Scientific",
        Maintainer = "Alan O'Callaghan <alan.ocallaghan@outlook.com>",
        RDataClass = "data.frame",
        DispatchClass = "Rds",
        stringsAsFactors = FALSE
    ),
    row.names = FALSE)
