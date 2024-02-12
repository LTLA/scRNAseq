
wrap_pmid = function(x)
 sprintf("<a href='https://pubmed.ncbi.nlm.nih.gov/%s'>%s</a>", x, x)


unbrack = function (x) 
gsub("\\{|\\}", "", x)

#' transform a bibtex file to data.frame
#' @import bibtex
#' @import dplyr
#' @param bibfile character(1) path to bibtex
#' @param has.pmids logical(1) if TRUE, href tags will be generated using the pmid attribute
#' of entries
#' @return data.frame with rownames given by the bibtex tags
#' @examples
#' bibref = system.file("scripts", "ref.bib", package="scRNAseq")
#' dem = bib_to_data.frame(bibref)
#' head(dem)
#' @export
bib_to_data.frame = function(bibfile, has.pmids=TRUE) {
 bibdat = bibtex::read.bib(bibfile)
 titles = sapply(bibdat, function(x) unbrack(x$title))
 if (has.pmids) {
   pmids = sapply(bibdat, function(x)
      wrap_pmid(x$pmid))
   return(data.frame(title=titles, pmid=pmids))
 }
 data.frame(title=titles)
}

#datf = bib_to_data.frame("ref.bib")

#' produce enriched data frame with pmid hrefs for studies in scRNAseq
#' @export
enhanceListDatasets = function() {
  bibref = system.file("scripts", "ref.bib", package="scRNAseq")
  bibdf = bib_to_data.frame("ref.bib")
  bibdf$tag = rownames(bibdf)
  ndf = as.data.frame(scRNAseq::listDatasets())
  ndf$tag = sub("@", "", ndf$Reference)
  nndf = dplyr::left_join(ndf, bibdf, by="tag")
  taxmap = c("8355"="Xenopus laevis", "9606"="Homo sapiens", "10090"="Mus musculus")
  data.frame(org=taxmap[as.character(nndf$Taxonomy)],
     part = nndf$Part, title=nndf$title, count=nndf$Number, pmid=nndf$pmid, call=nndf$Call)
}

#' use DT::datatable to present information in listDatasets() with linked HTML for pubmed entries
#' @import DT
#' @export
datasetsPage = function() {
  newdf = enhanceListDatasets()
  DT::datatable(newdf, escape=FALSE)
}
