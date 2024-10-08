% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleGeneAnnotation.R
\name{readGmt}
\alias{readGmt}
\title{This function reads in gene sets in .gmt format}
\usage{
readGmt(gmtFile, addDescr = FALSE)
}
\arguments{
\item{gmtFile}{\itemize{
\item a .gmt file containing gene sets, e.g., Hallmark of GO
}}

\item{addDescr}{\itemize{
\item include gene set description (2nd column in .gmt file) in
gene set name
}}
}
\value{
\itemize{
\item A named list of gene sets
}
}
\description{
This function reads in gene sets in .gmt format
}
