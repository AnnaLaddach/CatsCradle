% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cellTypesPerCellTypeMatrix}
\alias{cellTypesPerCellTypeMatrix}
\title{cellTypesPerCellTypeMatrix}
\format{
A matrix whose rows and columns correspond to
cell types.

\describe{
Each row of this matrix corresponds to a cell type.  On
that row we see the proportions of all neighbourhoods
surrounding cells of that cell type as regards the cell types
they contain.  In particular, each row sums to 1.
}

A matrix whose rows and columns correspond to
cell types.

\describe{
Each row of this matrix corresponds to a cell type.  On
that row we see the proportions of all neighbourhoods
surrounding cells of that cell type as regards the cell types
they contain.  In particular, each row sums to 1.
}
}
\source{
This is created from NBHDByCTMatrix and the clusters
using cellTypesPerCellTypeMatrix()

This is created from NBHDByCTMatrixExtended and the clusters
using cellTypesPerCellTypeMatrix()
}
\usage{
cellTypesPerCellTypeMatrix

cellTypesPerCellTypeMatrix
}
\description{
For each cell type, this matrix shows the fraction
of the neighbourhoods of that cell type composed of
each cell type.

For each cell type, this matrix shows the fraction
of the extended neighbourhoods of that cell type composed of
each cell type.
}
\keyword{datasets}
