% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleEdgeQualityControl.R
\name{edgeCutoffsByZScore}
\alias{edgeCutoffsByZScore}
\title{This finds edge cutoffs by z-score}
\usage{
edgeCutoffsByZScore(annEdges, zCutoff)
}
\arguments{
\item{annEdges}{\itemize{
\item a data frame with columns nodeA, nodeB, length
and cellTypePair as produced by edgeLengthsAndCellTypePairs.
}}

\item{zCutoff}{\itemize{
\item a numeric
}}
}
\value{
This returns a data frame with columns cellTypePair and
cutoff.
}
\description{
This finds edge cutoffs by z-score
}
\examples{
getExample = make.getExample()
centroids = getExample('centroids')
clusters = getExample('clusters')
delaunayNeighbours = getExample('delaunayNeighbours') 
annEdges =
    edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
cutoffDF = edgeCutoffsByZScore(annEdges,zCutoff=1.5)
}
