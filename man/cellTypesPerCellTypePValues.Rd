% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cellTypesPerCellTypePValues}
\alias{cellTypesPerCellTypePValues}
\title{cellTypesPerCellTypePValues
A symmetric matrix containing P values describing whether cell types are more
frequently neighbours than expected by chance.}
\format{
A matrix whose rows and columns correspond to
cell types.

\describe{
Rows and columns of this matrix correspond to a cell types. Matrix give p
values describing whether cell types are more frequently neighbours than
expected by chance.
}
}
\source{
This is created from delaunayNeighbours and the clusters
using computeNeighbourEnrichment()
}
\usage{
cellTypesPerCellTypePValues
}
\description{
cellTypesPerCellTypePValues
A symmetric matrix containing P values describing whether cell types are more
frequently neighbours than expected by chance.
}
\keyword{datasets}
