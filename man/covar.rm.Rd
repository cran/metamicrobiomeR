% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covar.rm.R
\docType{data}
\name{covar.rm}
\alias{covar.rm}
\title{Covariate data.}
\format{
A dataframe with 996 rows and 32 variables.
}
\source{
\href{https://gordonlab.wustl.edu/supplemental-data/supplemental-data-portal/subramanian-et-al-2014/}{Gordon Lab}
}
\usage{
data(covar.rm)
}
\description{
Monthly longitudinal clinical data of 50 infants from
birth to 2 years of life.
}
\examples{
data(covar.rm)
# Load KEGG pathway data
data(kegg.12)
# Comparison of pathway relative abundances for some first pathways of level 1 only
# and assuming crosssectional data (to save running time)
path1<-pathway.compare(pathtab=list(kegg.12[[1]][, 1:2]),
mapfile=covar.rm,sampleid="sampleid",pathsum="rel", stat.med="gamlss",
comvar="gender",adjustvar=c("age.sample","bf"), longitudinal="no",
p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)
taxcomtab.show(taxcomtab=path1$l1, sumvar="path",tax.lev="l2",
tax.select="none", showvar="genderMale", p.adjust.method="fdr",p.cutoff=1)
}
\references{
Subramanian et al. Nature. 2014 Jun 19; 510(7505): 417–421.
(\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4189846/}{PubMed})
}
\keyword{datasets}
