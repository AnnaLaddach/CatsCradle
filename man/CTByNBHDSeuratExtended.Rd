% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CTByNBHDSeuratExtended}
\alias{CTByNBHDSeuratExtended}
\title{CTByNBHDSeuratExtended}
\format{
A Seurat object consisting of 16 samples (cell types) and
4261 features (the neighbourhoods).

\describe{
This is a Seurat object created by taking t(NBHDByCTMatrixExtended)
as the counts.
}
}
\source{
Created from t(NBHDByCTMatrixExtended) by
computeNBHDVsCTSeurat()
}
\usage{
CTByNBHDSeuratExtended
}
\description{
A Seurat object computed from the transpose of NBHDByCTMatrixExtended.
Think of cell types "expressing" (being found in) neighbourhoods.
}
\keyword{datasets}
