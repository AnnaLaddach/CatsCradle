## ----setup, include = FALSE, warning = FALSE----------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    fig.dim = c(6,6),
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(Seurat)
library(CatsCradle)
ImageDimPlot(smallXenium, cols = "polychrome", size = 1)

## -----------------------------------------------------------------------------
centroids = GetTissueCoordinates(smallXenium)
delaunayNeighbours = computeNeighboursDelaunay(centroids)
head(delaunayNeighbours)

## -----------------------------------------------------------------------------
idx = (delaunayNeighbours$nodeA == 16307 |
       delaunayNeighbours$nodeB == 16307)
nbhd = unique(c(delaunayNeighbours$nodeA[idx],
                delaunayNeighbours$nodeB[idx]))
nbhd		

## -----------------------------------------------------------------------------
extendedNeighboursList = getExtendedNBHDs(delaunayNeighbours, 4)
extendedNeighbours = collapseExtendedNBHDs(extendedNeighboursList, 4)

## -----------------------------------------------------------------------------
idx = (extendedNeighbours$nodeA == 16307 |
       extendedNeighbours$nodeB == 16307)
nbhd = unique(c(extendedNeighbours$nodeA[idx],
                extendedNeighbours$nodeB[idx]))
length(nbhd)		

## -----------------------------------------------------------------------------
euclideanNeighbours = computeNeighboursEuclidean(centroids,
threshold=20)

