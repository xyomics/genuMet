% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simdata}
\alias{simdata}
\title{Simulated metabolomics data}
\format{Row headers:
\describe{
\item{date extracted}{The date the corresponding sample is processed.}
\item{Column}{The id of column that the sample is injected into. }
\item{injection order}{The injection order of the sample.}
}

Column headers:
\describe{
\item{metabolite}{Metabolite id.}
\item{method}{Experiment method.}
\item{HMDB id}{The id of the metabolite in HMDB.}
\item{mz}{Mass to charge ratio.}
\item{rt}{Retention time.}
\item{comment}{Comments for the metabolite.}
\item{processing compound id}{Processing compound id.}
}}
\usage{
data(simdata)
}
\description{
Simulated metabolomics dataset. The csv file consists of four parts.
The top-left of the table is blank. The top-right is row headers.
The down-left is column headers. And the down-right is metbolites signal.
}
\keyword{datasets}

