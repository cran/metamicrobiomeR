% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabsex4.R
\docType{data}
\name{tabsex4}
\alias{tabsex4}
\title{Combined data for meta-analysis.}
\format{
A dataframe with 701 rows and 23 variables.
}
\source{
\href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
}
\usage{
data(tabsex4)
}
\description{
Result outputs of differential abundance analysis using GAMLSS_BEZI
from "taxa.compare" function combined from
4 studies for meta-analysis. The comparison was between gender
adjusted for age of infants at sample collection.
}
\examples{
# load saved GAMLSS-BEZI results of four studies
# for the comparison of bacterial taxa relative abundance between
# genders adjusted for breastfeeding and infant age at sample collection
data(tabsex4)
#select only taxonomies of a small phylum for meta-analysis example
# (to save running time)
tlm<-tabsex4$id[grep("k__bacteria.p__fusobacteria",tabsex4$id)]
# meta-analysis
metab.sex<-meta.taxa(taxcomdat=tabsex4[tabsex4$id \%in\% tlm,],
summary.measure="RR", pool.var="id", studylab="study",
backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
#show results by table and plot
#phylum table
metatab.show(metatab=metab.sex$random,com.pooled.tab=tabsex4[tabsex4$id \%in\% tlm,],
tax.lev="l2",showvar="genderMale",p.cutoff.type="p", p.cutoff=1,display="table")
}
\references{
Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
(\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4189846/}{PubMed})

Bender et al. Sci Transl Med. 2016 Jul 27; 8(349): 349ra100.
(\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5301310/}{PubMed})

Pannaraj et al. JAMA Pediatr. 2017;90095(7):647–54.
(\href{https://pubmed.ncbi.nlm.nih.gov/28492938/}{PubMed})

Thompson et al. Front Cell Infect Microbiol. 2015;5:3.
(\href{https://pubmed.ncbi.nlm.nih.gov/25705611/}{PubMed})
}
\keyword{datasets}
