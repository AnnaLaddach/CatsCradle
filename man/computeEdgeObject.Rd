% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleLigandReceptorAnalysis.R
\name{computeEdgeObject}
\alias{computeEdgeObject}
\title{This function takes interactionResults and creates a seurat object where
each point represents an edge between cells, and spatial coordinates are the
centroids of edges between cells. The "expression matrix" is the
binarised presence/absence of an interaction (ligand receptor pair) on an edge.}
\usage{
computeEdgeObject(
  ligandReceptorResults,
  centroids,
  npcs = 10,
  returnType = "Seurat"
)
}
\arguments{
\item{ligandReceptorResults}{\itemize{
\item as returned by performLigandReceptorResultsAnalysis()
}}

\item{centroids}{\itemize{
\item a dataframe containing centroids
where rownames are cellnames and the first two columns
contain x and y coordinates respectively.
}}

\item{npcs}{\itemize{
\item number of pcs used for PCA, defaults to 10
}}

\item{returnType}{Determines whether to return a Seurat object or a
SpatialExperiment.  Will do the later if this is set to either SCE,
SingleCellExperiment or lower case versions of either.}
}
\value{
This returns a seurat object where
each point represents an edge between cells, and spatial coordinates are the
centroids of edges between cells. The "expression matrix" is the
binarised presence/absence of an interaction (ligand receptor pair) on an edge.
Depending on the parameter returnType, this can alternatively be returned as
a SpatialExperiment.
}
\description{
This function takes interactionResults and creates a seurat object where
each point represents an edge between cells, and spatial coordinates are the
centroids of edges between cells. The "expression matrix" is the
binarised presence/absence of an interaction (ligand receptor pair) on an edge.
}
\examples{
getExample = make.getExample()
centroids = getExample('centroids')
ligandReceptorResults = getExample('ligandReceptorResults')
edgeSeurat = computeEdgeObject(ligandReceptorResults, centroids)
}
