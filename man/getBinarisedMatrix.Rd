% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleLigandReceptorAnalysis.R
\name{getBinarisedMatrix}
\alias{getBinarisedMatrix}
\title{This functions retrieves an expression matrix from a
seurat object or SingleCellExperiment and binarises it.}
\usage{
getBinarisedMatrix(obj, cutoff = 0, layer = "count")
}
\arguments{
\item{obj}{\itemize{
\item a Seurat object or SingleCellExperiment to be
turned into a Seurat object
}}

\item{cutoff}{\itemize{
\item a cutoff for binarisation. Defaults to 0.
}}

\item{layer}{\itemize{
\item layer to fetch data from. Defaults to count.
}}
}
\value{
A binarised expression matrix where rows are genes and columns are
cells.
}
\description{
This functions retrieves an expression matrix from a
seurat object or SingleCellExperiment and binarises it.
}
