% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleNeighbourhoodAnalysis.R
\name{computeCellTypesPerCellTypeMatrix}
\alias{computeCellTypesPerCellTypeMatrix}
\title{For each cell type, this function looks at the neighbourhoods
around cells of that type and discovers the fractions of those
cells of each type.}
\usage{
computeCellTypesPerCellTypeMatrix(nbhdByCellType, cellTypes)
}
\arguments{
\item{nbhdByCellType}{\itemize{
\item A matrix whose rows are neighbourhoods
each denoted by the cell at their center, whose columns are
cell types, and whose entries are counts.
}}

\item{cellTypes}{\itemize{
\item named vector of cell types where names are each cell and
cell types are a factor
}}
}
\value{
A square matrix whose rownames and colnames are the
seurat_clusters as character strings.  Each row corresponds
to neighbourhoods around all cells of that type and the entries
give the fractions of those neighbourhoods occupied by cells
of each type.
}
\description{
For each cell type, this function looks at the neighbourhoods
around cells of that type and discovers the fractions of those
cells of each type.
}
\examples{
getExample = make.getExample()
NBHDByCTMatrix = getExample('NBHDByCTMatrix')
clusters = getExample('clusters')
cellTypesPerCellType = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,clusters)
}
