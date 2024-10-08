% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleNeighbourhoodAnalysis.R
\name{computeNBHDByCTMatrix}
\alias{computeNBHDByCTMatrix}
\title{This function computes a matrix where neighbourhoods are rows and
cell types are columns. The values in the matrix indicate the
number of cells of a given type within a neighbourhood.}
\usage{
computeNBHDByCTMatrix(spatialGraph, cellTypes)
}
\arguments{
\item{spatialGraph}{\itemize{
\item a spatial graph in neighbour list format.
}}

\item{cellTypes}{\itemize{
\item named vector of cell types where names are each cell and
cell types are a factor
}}
}
\value{
a matrix of neighbourhoods by cell types
}
\description{
This function computes a matrix where neighbourhoods are rows and
cell types are columns. The values in the matrix indicate the
number of cells of a given type within a neighbourhood.
}
\examples{
getExample = make.getExample()
clusters = getExample('clusters')
delaunayNeighbours = getExample('delaunayNeighbours')
NBHDByCTMatrix = computeNBHDByCTMatrix(delaunayNeighbours,clusters)
}
