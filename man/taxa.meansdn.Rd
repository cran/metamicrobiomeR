% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa.meansdn.R
\name{taxa.meansdn}
\alias{taxa.meansdn}
\title{Summarize abundance by group}
\usage{
taxa.meansdn(
  taxtab,
  sumvar,
  groupvar,
  percent.filter = 0.05,
  relabund.filter = 5e-05,
  othervar = "none"
)
}
\arguments{
\item{taxtab}{taxa/pathway  abundance table from phylum to species or any preferred highest taxa level.}

\item{sumvar}{main variable for summary}

\item{groupvar}{variable to be stratified.}

\item{percent.filter}{prevalence threshold (the percentage of number of samples the taxa/pathway available). Default is 0.05.}

\item{relabund.filter}{relative abundance threshold (the minimum of the average relative abundance for a taxa/pathway to be retained). Default is 0.00005.}

\item{othervar}{vector of variables that are not abundance variables to be summarized. Default is "none".}
}
\value{
table of mean, sd, count by group.
}
\description{
This function summarizes taxa/pathway abundance tables to provide mean, sd, count by groups.
}
\examples{
#Load summary tables of bacterial taxa relative abundance from Bangladesh data
data(taxtab6)
taxa.meansdn.rm<-taxa.meansdn(taxtab=taxtab6,sumvar="bf",groupvar="age.sample")
}
\keyword{abundance}
\keyword{summary}
