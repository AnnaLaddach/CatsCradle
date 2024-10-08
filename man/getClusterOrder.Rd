% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradleGeneAnnotation.R
\name{getClusterOrder}
\alias{getClusterOrder}
\title{This gets the clusters in their cannonical order}
\usage{
getClusterOrder(f)
}
\arguments{
\item{f}{\itemize{
\item a Seurat object with meta.data column seurat_clusters
or SingleCellExperiment to be turned into a Seurat object
}}
}
\value{
A vector of these unique values in order
}
\description{
This deals with skullduggery in which seurat_clusters
has been converted from a factor to a character or a
numeric.
}
\examples{
STranspose = make.getExample()('STranspose',toy=TRUE)
geneClusters = getClusterOrder(STranspose)
}
