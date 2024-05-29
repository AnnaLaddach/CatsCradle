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

## -----------------------------------------------------------------------------
library(ggplot2,quietly=TRUE)
h = 'HALLMARK_OXIDATIVE_PHOSPHORYLATION'
umap = FetchData(STranspose,c('umap_1','umap_2'))
idx = colnames(STranspose) %in% hallmark[[h]]
g = DimPlot(STranspose,cols='polychrome') +
    geom_point(data=umap[idx,],aes(x=umap_1,y=umap_2),color='black',size=2.7) +
    geom_point(data=umap[idx,],aes(x=umap_1,y=umap_2),color='green') +
    ggtitle(paste(h,'on gene clusters'))
print(g)
pValue = getSeuratSubsetClusteringPValue(STranspose,idx)
pValue

## -----------------------------------------------------------------------------
ImageDimPlot(smallXenium,cols='polychrome')

## -----------------------------------------------------------------------------
head(delaunayNeighbours,10)

## -----------------------------------------------------------------------------
cellTypesPerCellTypeMatrixExtended = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrixExtended,clusters)

cellTypesPerCellTypeGraphFromCellMatrix(cellTypesPerCellTypeMatrixExtended, minWeight = 0.05, colours = colours)

## -----------------------------------------------------------------------------
smallXenium$NBHDClusterExtended= 
  NBHDByCTSeuratExtended@active.ident
ImageDimPlot(smallXenium, group.by = c("NBHDClusterExtended"), 
             size = 1, cols = "polychrome")

## -----------------------------------------------------------------------------
agg = aggregateSeuratGeneExpression(smallXenium,extendedNeighbours,
                                    verbose=FALSE)
smallXenium$aggregateNBHDClusters = agg@active.ident
ImageDimPlot(smallXenium,group.by='aggregateNBHDClusters',cols='polychrome')

## -----------------------------------------------------------------------------
library(fossil,quietly=TRUE)
adjustedRandIndex = adj.rand.index(agg@active.ident,
                                   NBHDByCTSeuratExtended@active.ident)
adjustedRandIndex

## -----------------------------------------------------------------------------
ligandReceptorResults$interactionsOnEdges[1:10,1:10]

