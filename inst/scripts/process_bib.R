wrap_pmid = function(x)
 sprintf("<a href='https://pubmed.ncbi.nlm.nih.gov/%s'>%s</a>", x, x)


unbrack = function (x) 
gsub("\\{|\\}", "", x)

bib_to_datatable = function(bibfile) {
 bibdat = read.bib(bibfile)
 titles = sapply(bibdat, function(x) unbrack(x$title))
 pmids = sapply(bibdat, function(x)
  wrap_pmid(x$pmid))
 datf = data.frame(title=titles, pmid=pmids)
 datatable(datf, escape=FALSE)
}
