% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleLigandReceptorAnalysis.R
\name{makeSummedLRInteractionHeatmap}
\alias{makeSummedLRInteractionHeatmap}
\title{This function takes ligandReceptorResults and plots a heatmap of the total
number of ligand receptor interactions between clusters.}
\usage{
makeSummedLRInteractionHeatmap(
  ligandReceptorResults,
  clusters,
  type,
  logScale = TRUE
)
}
\arguments{
\item{ligandReceptorResults}{\itemize{
\item as returned by performLigandReceptorAnalysis()
}}

\item{clusters}{\itemize{
\item named vector of cell types where names are each cell and
clusters are a factor
}}

\item{type}{\itemize{
\item "total" or "mean" to plot raw total interactions or mean
interactions per edge.
}}

\item{logScale}{\itemize{
\item plot heatmap using log scale (defaults to TRUE)
}}
}
\value{
matrix of total ligand receptor interactions that underlies t
he heatmap.
}
\description{
This function takes ligandReceptorResults and plots a heatmap of the total
number of ligand receptor interactions between clusters.
}
\examples{
getExample = make.getExample()
clusters = getExample('clusters')
ligandReceptorResults = getExample('ligandReceptorResults')
cellTypePerCellTypeLigRecMatrix = 
makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "mean")
}
