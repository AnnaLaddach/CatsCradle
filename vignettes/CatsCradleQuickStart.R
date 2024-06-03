## ----setup, include = FALSE, warning = FALSE----------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    fig.dim = c(6,6),
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(Seurat,quietly=TRUE)
library(CatsCradle,quietly=TRUE)
DimPlot(S,cols='polychrome')

## -----------------------------------------------------------------------------
DimPlot(STranspose,cols='polychrome')

