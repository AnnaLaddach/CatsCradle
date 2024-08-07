% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CatsCradle.R
\name{transposeSeuratObject}
\alias{transposeSeuratObject}
\title{Create the transpose of a Seurat object}
\usage{
transposeSeuratObject(
  f,
  active.assay = "RNA",
  npcs = 30,
  dims = seq_len(20),
  res = 1,
  returnType = "Seurat"
)
}
\arguments{
\item{f}{\itemize{
\item a Seurat object
}}

\item{active.assay}{\itemize{
\item the assay to use.  Defaults to 'RNA'
}}

\item{npcs}{\itemize{
\item number of principal components, defaults to 30
}}

\item{dims}{\itemize{
\item dimensions to use for umap and nearest neighbors,
defaults to 1:20
}}

\item{res}{\itemize{
\item the clustering resolution, defaults to 1
}}

\item{returnType}{\itemize{
\item Will return a SingleCellExperiment if this is either
of SCE, SingleCellExperiment or their lower-case equivalents.  Otherwise,
returns a Seurat object
}}
}
\value{
A Seurat object or SingleCellExperiment
}
\description{
This takes a Seurat object f and creates a new Seurat object whose
expression matrix is the transpose of that of f.  This can
also be a SingleCellExperiment which will be converted to a
Seurat object
}
\examples{
STranspose = transposeSeuratObject(S)
STransposeAsSCE = transposeSeuratObject(S,returnType='SCE')
}
